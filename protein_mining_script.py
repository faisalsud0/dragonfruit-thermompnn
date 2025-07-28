#!/usr/bin/env python3
"""
Plant Thermostable Stress Response Protein Mining Script

This script searches NCBI databases for heat shock proteins and ROS detoxifying proteins
across a list of plant species, collecting sequence and structural information suitable
for AI models like AlphaFold.

Requirements:
- biopython
- pandas
- requests
- time
- csv

Install with: pip install biopython pandas requests
"""

import pandas as pd
import csv
import time
import requests
import json
import re
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from typing import List, Dict, Tuple, Optional
import logging
from datetime import datetime
import os

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('protein_mining.log'),
        logging.StreamHandler()
    ]
)

class ProteinMiner:
    def __init__(self, email: str, api_key: str = None):
        """
        Initialize the protein mining system.
        
        Args:
            email: Your email for NCBI API access (required)
            api_key: NCBI API key for higher rate limits (optional but recommended)
        """
        # Set up NCBI Entrez
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        
        # Define protein search terms
        self.heat_shock_terms = [
            "heat shock protein", "HSP", "heat stress protein", "molecular chaperone",
            "chaperone", "DnaK", "DnaJ", "GroEL", "GroES", "ClpX", "ClpP",
            "small heat shock protein", "sHSP", "HSP70", "HSP60", "HSP90", "HSP100",
            "HSP20", "HSP27", "HSP40", "chaperonin"
        ]
        
        self.ros_detox_terms = [
            "superoxide dismutase", "SOD", "catalase", "CAT", "peroxidase", "POD",
            "glutathione peroxidase", "GPX", "ascorbate peroxidase", "APX",
            "thioredoxin", "TRX", "glutaredoxin", "GRX", "peroxiredoxin", "PRX",
            "glutathione reductase", "GR", "dehydroascorbate reductase", "DHAR",
            "monodehydroascorbate reductase", "MDHAR", "glutathione S-transferase", "GST",
            "aldehyde dehydrogenase", "ALDH", "glyoxalase", "GLY", "cytochrome c peroxidase"
        ]
        
        # Rate limiting
        self.request_delay = 0.34  # ~3 requests per second (NCBI limit without API key)
        if api_key:
            self.request_delay = 0.1  # 10 requests per second with API key
    
    def load_species_data(self, csv_file: str) -> pd.DataFrame:
        """Load plant species data from CSV file."""
        try:
            # Read CSV, skipping the first empty row
            df = pd.read_csv(csv_file, skiprows=1)
            
            # Print column names for debugging
            logging.info(f"CSV columns found: {list(df.columns)}")
            
            # Find the species column (it might have different name)
            species_col = None
            for col in df.columns:
                if col and ('species' in col.lower() or 'scientific' in col.lower()):
                    species_col = col
                    break
            
            if species_col is None:
                # If no species column found, assume it's the 5th column (index 4) based on your CSV structure
                if len(df.columns) >= 5:
                    species_col = df.columns[4]  # 'species' column
                else:
                    raise ValueError("Could not find species column in CSV")
            
            logging.info(f"Using column '{species_col}' for species names")
            
            # Clean and validate species names
            df = df.dropna(subset=[species_col])
            df[species_col] = df[species_col].astype(str).str.strip()
            
            # Remove entries with generic terms
            df = df[~df[species_col].str.contains(r'sp\.|species|var\.|cf\.', case=False, na=False)]
            
            # Rename the species column for consistency
            df = df.rename(columns={species_col: 'species'})
            
            logging.info(f"Loaded {len(df)} valid species entries")
            return df
        except Exception as e:
            logging.error(f"Error loading CSV file: {e}")
            raise
    
    def create_search_query(self, species: str, protein_terms: List[str]) -> str:
        """Create NCBI search query for specific species and protein types."""
        # Try multiple query strategies
        species_clean = species.replace(" ", "[Organism] AND ")
        
        # Create a broader search - combine terms with OR, but make it less restrictive
        protein_query = " OR ".join([f'{term}' for term in protein_terms])
        
        # Try simpler query first
        query = f'({species_clean}[Organism]) AND ({protein_query})'
        return query
    
    def search_ncbi_proteins(self, species: str, protein_type: str) -> List[Dict]:
        """
        Search NCBI Protein database for specific protein types in a species.
        
        Args:
            species: Scientific name of the species
            protein_type: 'heat_shock' or 'ros_detox'
        
        Returns:
            List of protein records with metadata
        """
        protein_records = []
        
        try:
            # Select appropriate search terms
            if protein_type == 'heat_shock':
                terms = self.heat_shock_terms[:5]  # Use fewer, more common terms
            elif protein_type == 'ros_detox':
                terms = self.ros_detox_terms[:5]  # Use fewer, more common terms
            else:
                raise ValueError("protein_type must be 'heat_shock' or 'ros_detox'")
            
            # Create search query
            query = self.create_search_query(species, terms)
            logging.info(f"Searching {species} for {protein_type} proteins with query: {query}")
            
            # Search NCBI
            search_handle = Entrez.esearch(
                db="protein",
                term=query,
                retmax=100,  # Reduced limit
                sort="relevance"
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            protein_ids = search_results["IdList"]
            
            if not protein_ids:
                # Try a simpler search with just the species and basic terms
                simple_terms = ['heat shock', 'HSP'] if protein_type == 'heat_shock' else ['superoxide dismutase', 'catalase', 'peroxidase']
                simple_query = f'{species}[Organism] AND ({" OR ".join(simple_terms)})'
                logging.info(f"Trying simpler query: {simple_query}")
                
                search_handle = Entrez.esearch(
                    db="protein",
                    term=simple_query,
                    retmax=50,
                    sort="relevance"
                )
                search_results = Entrez.read(search_handle)
                search_handle.close()
                protein_ids = search_results["IdList"]
            
            if not protein_ids:
                logging.info(f"No {protein_type} proteins found for {species}")
                return protein_records
            
            logging.info(f"Found {len(protein_ids)} {protein_type} proteins for {species}")
            
            # Fetch detailed information in batches
            batch_size = 20
            for i in range(0, len(protein_ids), batch_size):
                batch_ids = protein_ids[i:i + batch_size]
                
                # Rate limiting
                time.sleep(self.request_delay)
                
                try:
                    # Fetch protein records
                    fetch_handle = Entrez.efetch(
                        db="protein",
                        id=batch_ids,
                        rettype="gb",
                        retmode="xml"
                    )
                    records = Entrez.read(fetch_handle)
                    fetch_handle.close()
                    
                    # Process each record
                    for record in records:
                        protein_info = self.extract_protein_info(record, species, protein_type)
                        if protein_info:
                            protein_records.append(protein_info)
                
                except Exception as e:
                    logging.warning(f"Error fetching batch for {species}: {e}")
                    continue
        
        except Exception as e:
            logging.error(f"Error searching {species} for {protein_type}: {e}")
        
        return protein_records
    
    def extract_protein_info(self, record: Dict, species: str, protein_type: str) -> Optional[Dict]:
        """Extract relevant information from NCBI protein record."""
        try:
            # Get basic information
            protein_info = {
                'species': species,
                'protein_type': protein_type,
                'accession': record.get('GBSeq_accession-version', ''),
                'gi': record.get('GBSeq_primary-accession', ''),
                'length': int(record.get('GBSeq_length', 0)),
                'definition': record.get('GBSeq_definition', ''),
                'organism': record.get('GBSeq_organism', ''),
                'sequence': record.get('GBSeq_sequence', '').upper(),
                'create_date': record.get('GBSeq_create-date', ''),
                'update_date': record.get('GBSeq_update-date', ''),
                'topology': record.get('GBSeq_topology', ''),
                'division': record.get('GBSeq_division', ''),
            }
            
            # Extract features and annotations
            features = record.get('GBSeq_feature-table', [])
            protein_info['features'] = []
            
            for feature in features:
                feature_type = feature.get('GBFeature_key', '')
                if feature_type in ['CDS', 'Protein', 'mat_peptide', 'sig_peptide']:
                    feature_info = {
                        'type': feature_type,
                        'location': feature.get('GBFeature_location', ''),
                        'qualifiers': {}
                    }
                    
                    # Extract qualifiers
                    qualifiers = feature.get('GBFeature_quals', [])
                    for qual in qualifiers:
                        qual_name = qual.get('GBQualifier_name', '')
                        qual_value = qual.get('GBQualifier_value', '')
                        feature_info['qualifiers'][qual_name] = qual_value
                    
                    protein_info['features'].append(feature_info)
            
            # Extract references
            references = record.get('GBSeq_references', [])
            protein_info['references'] = []
            
            for ref in references[:3]:  # Limit to first 3 references
                ref_info = {
                    'title': ref.get('GBReference_title', ''),
                    'authors': ref.get('GBReference_authors', []),
                    'journal': ref.get('GBReference_journal', ''),
                    'pubmed': ref.get('GBReference_pubmed', '')
                }
                protein_info['references'].append(ref_info)
            
            # Calculate basic sequence properties
            if protein_info['sequence']:
                protein_info['molecular_weight'] = self.calculate_molecular_weight(protein_info['sequence'])
                protein_info['amino_acid_composition'] = self.get_aa_composition(protein_info['sequence'])
            
            return protein_info
        
        except Exception as e:
            logging.warning(f"Error extracting protein info: {e}")
            return None
    
    def calculate_molecular_weight(self, sequence: str) -> float:
        """Calculate approximate molecular weight of protein sequence."""
        # Approximate molecular weights of amino acids
        aa_weights = {
            'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
            'Q': 146.15, 'E': 147.13, 'G': 75.07, 'H': 155.16, 'I': 131.17,
            'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
            'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
        }
        
        weight = sum(aa_weights.get(aa, 0) for aa in sequence.upper())
        # Subtract water molecules for peptide bonds
        weight -= (len(sequence) - 1) * 18.015
        return round(weight, 2)
    
    def get_aa_composition(self, sequence: str) -> Dict[str, float]:
        """Calculate amino acid composition as percentages."""
        if not sequence:
            return {}
        
        counts = {}
        for aa in sequence.upper():
            counts[aa] = counts.get(aa, 0) + 1
        
        total = len(sequence)
        composition = {aa: round((count / total) * 100, 2) for aa, count in counts.items()}
        return composition
    
    def mine_species_proteins(self, species: str) -> List[Dict]:
        """Mine both heat shock and ROS detoxifying proteins for a species."""
        all_proteins = []
        
        # Search for heat shock proteins
        heat_shock_proteins = self.search_ncbi_proteins(species, 'heat_shock')
        all_proteins.extend(heat_shock_proteins)
        
        # Small delay between searches
        time.sleep(self.request_delay)
        
        # Search for ROS detoxifying proteins
        ros_detox_proteins = self.search_ncbi_proteins(species, 'ros_detox')
        all_proteins.extend(ros_detox_proteins)
        
        return all_proteins
    
    def save_results(self, results: List[Dict], output_file: str):
        """Save mining results to CSV file."""
        if not results:
            logging.warning("No results to save")
            return
        
        # Flatten the data for CSV export
        flattened_results = []
        
        for protein in results:
            base_info = {
                'species': protein['species'],
                'protein_type': protein['protein_type'],
                'accession': protein['accession'],
                'gi': protein['gi'],
                'length': protein['length'],
                'definition': protein['definition'],
                'organism': protein['organism'],
                'sequence': protein['sequence'],
                'create_date': protein['create_date'],
                'update_date': protein['update_date'],
                'molecular_weight': protein.get('molecular_weight', ''),
                'topology': protein['topology'],
                'division': protein['division']
            }
            
            # Add amino acid composition
            aa_comp = protein.get('amino_acid_composition', {})
            for aa, percentage in aa_comp.items():
                base_info[f'aa_{aa}_percent'] = percentage
            
            # Add feature information
            features = protein.get('features', [])
            base_info['num_features'] = len(features)
            
            # Add reference information
            references = protein.get('references', [])
            base_info['num_references'] = len(references)
            if references:
                base_info['first_reference_title'] = references[0].get('title', '')
                base_info['first_reference_journal'] = references[0].get('journal', '')
                base_info['first_reference_pubmed'] = references[0].get('pubmed', '')
            
            flattened_results.append(base_info)
        
        # Save to CSV
        df = pd.DataFrame(flattened_results)
        df.to_csv(output_file, index=False)
        logging.info(f"Results saved to {output_file}")
    
    def save_fasta_sequences(self, results: List[Dict], output_file: str):
        """Save protein sequences in FASTA format for AlphaFold/ML models."""
        with open(output_file, 'w') as f:
            for protein in results:
                if protein['sequence']:
                    header = f">{protein['accession']}|{protein['species']}|{protein['protein_type']}|{protein['definition'][:100]}"
                    f.write(header + '\n')
                    # Write sequence in 80-character lines
                    sequence = protein['sequence']
                    for i in range(0, len(sequence), 80):
                        f.write(sequence[i:i+80] + '\n')
        
        logging.info(f"FASTA sequences saved to {output_file}")
    
    def test_species_availability(self, species: str) -> Dict:
        """Test what proteins are available for a species in NCBI."""
        try:
            # Simple search for any proteins from this species
            simple_query = f'{species}[Organism]'
            
            search_handle = Entrez.esearch(
                db="protein",
                term=simple_query,
                retmax=10
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            total_proteins = int(search_results["Count"])
            protein_ids = search_results["IdList"]
            
            logging.info(f"Species {species}: {total_proteins} total proteins in NCBI")
            
            if protein_ids:
                # Get a sample protein to see what data looks like
                fetch_handle = Entrez.efetch(
                    db="protein",
                    id=protein_ids[0],
                    rettype="gb",
                    retmode="xml"
                )
                records = Entrez.read(fetch_handle)
                fetch_handle.close()
                
                if records:
                    sample_protein = records[0]
                    logging.info(f"Sample protein: {sample_protein.get('GBSeq_definition', 'No definition')}")
            
            return {
                'species': species,
                'total_proteins': total_proteins,
                'sample_ids': protein_ids[:5]
            }
        
        except Exception as e:
            logging.error(f"Error testing species {species}: {e}")
            return {'species': species, 'total_proteins': 0, 'sample_ids': []}

    def run_species_test(self, csv_file: str, num_species: int = 5):
        """Test a few species to see what's available."""
        species_df = self.load_species_data(csv_file)
        
        logging.info("=== TESTING SPECIES AVAILABILITY ===")
        
        for i in range(min(num_species, len(species_df))):
            species = species_df.iloc[i]['species']
            result = self.test_species_availability(species)
            time.sleep(self.request_delay)
        
        logging.info("=== END SPECIES TEST ===")

    def run_complete_mining(self, csv_file: str, output_prefix: str = "protein_mining_results", 
                           max_species: int = None, start_index: int = 0):
        """
        Run complete protein mining pipeline.
        
        Args:
            csv_file: Input CSV file with species data
            output_prefix: Prefix for output files
            max_species: Maximum number of species to process (None for all)
            start_index: Starting index for resuming interrupted runs
        """
        # Load species data
        species_df = self.load_species_data(csv_file)
        
        # Limit species if specified
        if max_species:
            species_df = species_df.iloc[start_index:start_index + max_species]
        else:
            species_df = species_df.iloc[start_index:]
        
        all_results = []
        processed_count = 0
        
        # Create progress file for resuming
        progress_file = f"{output_prefix}_progress.txt"
        
        for index, row in species_df.iterrows():
            species = row['species']
            
            try:
                logging.info(f"Processing species {processed_count + 1}/{len(species_df)}: {species}")
                
                # Mine proteins for this species
                species_proteins = self.mine_species_proteins(species)
                all_results.extend(species_proteins)
                
                processed_count += 1
                
                # Save progress
                with open(progress_file, 'w') as f:
                    f.write(f"Last processed: {species}\n")
                    f.write(f"Index: {start_index + processed_count}\n")
                    f.write(f"Total proteins found so far: {len(all_results)}\n")
                
                # Save intermediate results every 50 species
                if processed_count % 50 == 0:
                    intermediate_file = f"{output_prefix}_intermediate_{processed_count}.csv"
                    self.save_results(all_results, intermediate_file)
                    
                    intermediate_fasta = f"{output_prefix}_intermediate_{processed_count}.fasta"
                    self.save_fasta_sequences(all_results, intermediate_fasta)
                    
                    logging.info(f"Saved intermediate results after {processed_count} species")
                
                # Longer delay between species to be respectful to NCBI
                time.sleep(self.request_delay * 2)
                
            except Exception as e:
                logging.error(f"Error processing {species}: {e}")
                continue
        
        # Save final results
        final_csv = f"{output_prefix}_final.csv"
        final_fasta = f"{output_prefix}_final.fasta"
        
        self.save_results(all_results, final_csv)
        self.save_fasta_sequences(all_results, final_fasta)
        
        # Summary statistics
        logging.info(f"Mining complete!")
        logging.info(f"Processed {processed_count} species")
        logging.info(f"Found {len(all_results)} total proteins")
        
        # Count by protein type
        heat_shock_count = sum(1 for p in all_results if p['protein_type'] == 'heat_shock')
        ros_detox_count = sum(1 for p in all_results if p['protein_type'] == 'ros_detox')
        
        logging.info(f"Heat shock proteins: {heat_shock_count}")
        logging.info(f"ROS detoxifying proteins: {ros_detox_count}")
        
        return all_results


def main():
    """Main execution function with example usage."""
    
    # CONFIGURATION - MODIFY THESE VALUES
    EMAIL = "alrajihsaad@gmail.com"  # REQUIRED: Your email for NCBI
    API_KEY = "f8d15735b147912be28808a1d150f65c2207"  # OPTIONAL: Your NCBI API key for higher rate limits
    CSV_FILE = r"C:\Users\alraj\Downloads\Data Collection - Litereature Info.csv"  # Your input CSV file
    OUTPUT_PREFIX = "plant_stress_proteins"  # Prefix for output files
    
    # Initialize the miner
    miner = ProteinMiner(email=EMAIL, api_key=API_KEY)
    
    # Run the complete mining pipeline for ALL species
    print("Starting protein mining for all 1,073 species...")
    print("This will take several hours but saves progress every 50 species.")
    results = miner.run_complete_mining(
        csv_file=CSV_FILE,
        output_prefix=OUTPUT_PREFIX,
        max_species=None,  # Process all 1,073 species
        start_index=0     # Use this to resume from a specific point
    )
    
    print(f"Mining completed! Found {len(results)} proteins total.")
    print(f"Check the output files: {OUTPUT_PREFIX}_final.csv and {OUTPUT_PREFIX}_final.fasta")


if __name__ == "__main__":
    main()