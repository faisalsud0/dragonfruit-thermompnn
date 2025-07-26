# dragonfruit-thermompnn

## Installation Steps

### 1. Install Miniconda

1. Download Miniconda from: https://docs.conda.io/en/latest/miniconda.html
2. Choose "Miniconda3 Windows 64-bit" installer
3. Run the installer with these settings:
    - ✅ Check "Add Miniconda3 to PATH environment variable"
    - ✅ Check "Register Miniconda3 as my default Python"
4. Restart your computer after installation

### 2. Verify Conda Installation

Open the terminal and test:

```bash
conda --version

```

You should see something like `conda 25.x.x`

### 3. Create Environment

```bash
conda create -n dragonfruit python=3.14
conda activate dragonfruit
python --version

```

### 4. Clone the Repositories

```bash
git clone https://github.com/dauparas/ProteinMPNN.git
git clone https://github.com/Kuhlman-Lab/ThermoMPNN.git

```

### 5. Verify Downloads

You should now have:

```
dragonfruit-thermompnn/
├── ProteinMPNN/
├── ThermoMPNN/
└── README.md

```

### 6. Activate Environment

Activate your environment before working:

```bash
conda activate dragonfruit

```

### 7. Installing ProteinMPNN and ThermoMPNN Dependencies

**Install PyTorch with CUDA support:**

```bash
conda install pytorch pytorch-cuda=11.8 pytorch-lightning -c pytorch -c nvidia -c conda-forge

```

**Install all other dependencies:**

```bash
conda install numpy joblib omegaconf pandas tqdm mmseqs2 wandb biopython -c bioconda -c conda-forge -c anaconda

```

**Check GPU recognition (should return True)**

```bash
python -c 'import torch; print(torch.cuda.is_available())'
cd ..

```

### 8. Test Installation

From your root repo directory, run `python test_installation.py`