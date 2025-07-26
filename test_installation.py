#!/usr/bin/env python3
"""
Test Installation Script for dragonfruit-thermompnn
Tests all dependencies and components for ProteinMPNN and ThermoMPNN
"""

import sys
import os
import subprocess
from pathlib import Path

def print_header(text):
    """Print a formatted header"""
    print(f"\n{'='*60}")
    print(f" {text}")
    print(f"{'='*60}")

def print_test(test_name, status, details=""):
    """Print test results in a formatted way"""
    status_symbol = "✅" if status else "❌"
    print(f"{status_symbol} {test_name}")
    if details:
        print(f"   {details}")

def test_python_version():
    """Test Python version"""
    print_header("PYTHON VERSION CHECK")
    
    version = sys.version_info
    print(f"Python version: {version.major}.{version.minor}.{version.micro}")
    
    # Check if Python 3.10+
    if version.major == 3 and version.minor >= 10:
        print_test("Python Version", True, f"Python {version.major}.{version.minor}.{version.micro} ✓")
        return True
    else:
        print_test("Python Version", False, f"Need Python 3.10+, got {version.major}.{version.minor}.{version.micro}")
        return False

def test_conda_environment():
    """Test conda environment"""
    print_header("CONDA ENVIRONMENT CHECK")
    
    conda_env = os.environ.get('CONDA_DEFAULT_ENV', 'None')
    print_test("Conda Environment", conda_env == 'dragonfruit', f"Active environment: {conda_env}")
    
    return conda_env == 'dragonfruit'

def test_core_imports():
    """Test core Python package imports"""
    print_header("CORE PACKAGE IMPORTS")
    
    packages = [
        'numpy',
        'pandas', 
        'torch',
        'joblib',
        'omegaconf',
        'tqdm',
        'Bio',  # biopython
        'wandb'
    ]
    
    results = []
    for package in packages:
        try:
            __import__(package)
            print_test(f"Import {package}", True)
            results.append(True)
        except ImportError as e:
            print_test(f"Import {package}", False, f"Error: {e}")
            results.append(False)
    
    return all(results)

def test_pytorch_cuda():
    """Test PyTorch CUDA availability"""
    print_header("PYTORCH CUDA CHECK")
    
    try:
        import torch
        
        print(f"PyTorch version: {torch.__version__}")
        
        cuda_available = torch.cuda.is_available()
        print_test("CUDA Available", cuda_available)
        
        if cuda_available:
            device_count = torch.cuda.device_count()
            current_device = torch.cuda.current_device()
            device_name = torch.cuda.get_device_name(current_device)
            
            print(f"   GPU Count: {device_count}")
            print(f"   Current Device: {current_device}")
            print(f"   Device Name: {device_name}")
            
            # Test basic CUDA operations
            try:
                x = torch.tensor([1.0, 2.0]).cuda()
                y = x * 2
                print_test("CUDA Operations", True, "Basic tensor operations work")
                return True
            except Exception as e:
                print_test("CUDA Operations", False, f"Error: {e}")
                return False
        else:
            print("   No CUDA devices found - will use CPU")
            return False
            
    except ImportError:
        print_test("PyTorch Import", False, "PyTorch not installed")
        return False

def test_repository_structure():
    """Test if repositories are cloned correctly"""
    print_header("REPOSITORY STRUCTURE CHECK")
    
    current_dir = Path.cwd()
    print(f"Current directory: {current_dir}")
    
    # Check for ProteinMPNN
    proteinmpnn_path = current_dir / "ProteinMPNN"
    proteinmpnn_exists = proteinmpnn_path.exists() and proteinmpnn_path.is_dir()
    print_test("ProteinMPNN directory", proteinmpnn_exists, f"Path: {proteinmpnn_path}")
    
    # Check for ThermoMPNN  
    thermompnn_path = current_dir / "ThermoMPNN"
    thermompnn_exists = thermompnn_path.exists() and thermompnn_path.is_dir()
    print_test("ThermoMPNN directory", thermompnn_exists, f"Path: {thermompnn_path}")
    
    # Check for key files in ProteinMPNN
    if proteinmpnn_exists:
        key_files = [
            proteinmpnn_path / "protein_mpnn_run.py",
            proteinmpnn_path / "helper_scripts",
        ]
        
        for file_path in key_files:
            exists = file_path.exists()
            print_test(f"ProteinMPNN/{file_path.name}", exists)
    
    # Check for key files in ThermoMPNN
    if thermompnn_exists:
        key_files = [
            thermompnn_path / "thermompnn",
            thermompnn_path / "setup.py",
        ]
        
        for file_path in key_files:
            exists = file_path.exists()
            print_test(f"ThermoMPNN/{file_path.name}", exists)
    
    return proteinmpnn_exists and thermompnn_exists

def test_external_tools():
    """Test external command-line tools"""
    print_header("EXTERNAL TOOLS CHECK")
    
    tools = ['mmseqs', 'git']
    results = []
    
    for tool in tools:
        try:
            result = subprocess.run([tool, '--version'], 
                                  capture_output=True, 
                                  text=True, 
                                  timeout=10)
            if result.returncode == 0:
                version = result.stdout.split('\n')[0]
                print_test(f"{tool}", True, f"Version: {version}")
                results.append(True)
            else:
                print_test(f"{tool}", False, "Not found or error")
                results.append(False)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            print_test(f"{tool}", False, "Not found in PATH")
            results.append(False)
    
    return all(results)

def test_simple_torch_model():
    """Test a simple PyTorch model to ensure everything works"""
    print_header("PYTORCH MODEL TEST")
    
    try:
        import torch
        import torch.nn as nn
        
        # Create a simple model
        model = nn.Linear(10, 1)
        
        # Create sample data
        x = torch.randn(5, 10)
        
        # Forward pass
        with torch.no_grad():
            output = model(x)
        
        print_test("Simple Model Creation", True, f"Output shape: {output.shape}")
        
        # Test CUDA if available
        if torch.cuda.is_available():
            try:
                model_cuda = model.cuda()
                x_cuda = x.cuda()
                output_cuda = model_cuda(x_cuda)
                print_test("CUDA Model Test", True, "GPU computation successful")
                return True
            except Exception as e:
                print_test("CUDA Model Test", False, f"Error: {e}")
                return False
        else:
            print_test("CUDA Model Test", False, "CUDA not available - using CPU")
            return True
            
    except Exception as e:
        print_test("PyTorch Model Test", False, f"Error: {e}")
        return False

def run_all_tests():
    """Run all tests and provide summary"""
    print_header("DRAGONFRUIT-THERMOMPNN INSTALLATION TEST")
    print("Testing all components required for ProteinMPNN and ThermoMPNN...")
    
    tests = [
        ("Python Version", test_python_version),
        ("Conda Environment", test_conda_environment), 
        ("Core Imports", test_core_imports),
        ("PyTorch CUDA", test_pytorch_cuda),
        ("Repository Structure", test_repository_structure),
        ("External Tools", test_external_tools),
        ("PyTorch Model", test_simple_torch_model)
    ]
    
    results = {}
    
    for test_name, test_func in tests:
        try:
            results[test_name] = test_func()
        except Exception as e:
            print_test(f"{test_name} (CRASHED)", False, f"Unexpected error: {e}")
            results[test_name] = False
    
    # Summary
    print_header("TEST SUMMARY")
    
    passed = sum(results.values())
    total = len(results)
    
    for test_name, result in results.items():
        status = "PASS" if result else "FAIL"
        symbol = "✅" if result else "❌"
        print(f"{symbol} {test_name}: {status}")
    
    print(f"\nOverall: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 ALL TESTS PASSED! Your installation is ready for ProteinMPNN and ThermoMPNN!")
    else:
        print(f"\n⚠️  {total - passed} test(s) failed. Please fix the issues above before proceeding.")
        print("\nCommon fixes:")
        print("- Make sure you're in the 'dragonfruit' conda environment")
        print("- Install missing packages with conda install")
        print("- Clone missing repositories with git clone")
        print("- Ensure CUDA drivers are installed for GPU support")
    
    return passed == total

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)