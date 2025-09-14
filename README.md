# Basilisk-101nano: Hands-On CFD Simulations Using Basilisk C

[![Issues](https://img.shields.io/github/issues/comphy-lab/Basilisk-101)](https://github.com/comphy-lab/Basilisk-101/issues)
[![License](https://img.shields.io/github/license/comphy-lab/Basilisk-101)](https://github.com/comphy-lab/Basilisk-101/blob/main/LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/comphy-lab/Basilisk-101)](https://github.com/comphy-lab/Basilisk-101/commits/main)
[![Course](https://img.shields.io/badge/Course-Nano%20Version-blue)](https://comphy-lab.org/teaching/)
[![Basilisk](https://img.shields.io/badge/Basilisk-Compatible-green)](http://basilisk.fr/)
![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/comphy-lab/Basilisk-101?utm_source=oss&utm_medium=github&utm_campaign=comphy-lab%2FBasilisk-101&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)

**A compact, hands-on introduction to Basilisk C for computational fluid dynamics**

This is the **nano version** of the Basilisk-101 course - a fast-tracked, 3-hour intensive workshop focused on practical implementation rather than lectures. Perfect for researchers who want to quickly get up and running with adaptive mesh CFD simulations.

## Course Description

This **3-hour intensive workshop** provides a hands-on introduction to computational fluid dynamics using Basilisk C. Unlike traditional lecture-heavy courses, this nano version emphasizes **learning by doing** - you'll spend most of your time coding, running simulations, and exploring physics in real-time.

### What You'll Learn

- **Hour 1: Foundations** - Think before you compute! Basic transport equations and your first Basilisk simulation
- **Hour 2: Interface Dynamics** - Volume-of-Fluid methods and drop impact physics with live parameter variations
- **Hour 3: Coating Applications** - Contact line dynamics and Landau-Levich dip coating with immediate visual feedback

### Key Features

- **Hands-on focused**: 80% coding, 20% theory
- **Interactive exploration**: Modify parameters and see immediate results
- **Physics-first approach**: Understanding the physics before numerical implementation
- **Production-ready code**: Learn with actual research-grade simulations

### Course Format

- **Duration**: 3 hours intensive
- **Format**: Workshop-style with live coding
- **Prerequisites**:
  - Basic fluid mechanics knowledge
  - Some programming experience (any language)
  - Laptop with C compiler capability
- **Materials**: All code examples and exercises included in this repository

Perfect for researchers, graduate students, and engineers who want to quickly start using adaptive mesh CFD for their own projects.

## Repository Structure

```
├── basilisk/src/: Core Basilisk CFD library (reference only)
├── postProcess/: Post-processing tools for visualization
├── src-local/: Custom header files extending Basilisk functionality
├── testCases/: Workshop exercises and examples
    ├── 1-*: Foundations - Transport equations & basic CFD
    ├── 2-*: Flow dynamics - Navier-Stokes applications
    ├── 3-DropImpactOnSolids.c: Interface dynamics workshop
    ├── 4-DipCoating-Withdrawal.c: Landau-Levich coating physics
    └── 4-DipCoating-Plunging.c: Contact line dynamics
```

### Key Workshop Files

The nano course focuses on these essential examples:

- **Hour 1**: `1-conduction-*.c` - Foundation physics and Basilisk syntax
- **Hour 2**: `3-DropImpactOnSolids.c` - Volume-of-Fluid and adaptive mesh refinement
- **Hour 3**: `4-DipCoating-*.c` - Industrial coating applications with contact line physics

Each file is fully commented and designed for live parameter exploration during the workshop.

## Quick Start (For Workshop Participants)

If you're attending a nano workshop, follow these streamlined steps:

1. **Clone this repository**:
   ```bash
   git clone https://github.com/comphy-lab/Basilisk-101nano.git
   cd Basilisk-101nano
   ```

2. **Quick Basilisk setup**:
   ```bash
   ./reset_install_requirements.sh
   source .project_config
   ```

3. **Test your installation**:
   ```bash
   cd testCases
   CFLAGS=-DDISPLAY=-1 make 1-conduction-simple.tst
   ```

4. **Ready to code!** Your instructor will guide you through the examples.

For detailed installation instructions or troubleshooting, see the full installation section below.

## Installation

Basilisk requires a C99-compliant compiler and GNU make. Installation can be done in two ways:

### Global Installation

#### Using darcs (recommended)
```bash
sudo apt install darcs make gawk
darcs clone http://basilisk.fr/basilisk
cd basilisk/src

# For Linux/Ubuntu users (preferred mode of operation)
ln -s config.gcc config

# For Mac users
# ln -s config.osx config

make
```

#### Using a tarball
```bash
wget http://basilisk.fr/basilisk/basilisk.tar.gz
tar xzf basilisk.tar.gz
cd basilisk/src

# For Linux/Ubuntu users (preferred mode of operation)
ln -s config.gcc config

# For Mac users
# ln -s config.osx config

make
```

#### Add to your shell configuration (.bashrc or .zshrc)
```bash
echo "export BASILISK=$PWD" >> ~/.bashrc
echo 'export PATH=$PATH:$BASILISK' >> ~/.bashrc
```
Or for zsh users:
```bash
echo "export BASILISK=$PWD" >> ~/.zshrc
echo 'export PATH=$PATH:$BASILISK' >> ~/.zshrc
```

### Repository Level Installation

For project-specific installations, you can use the provided `reset_install_requirements.sh` script which:
- Installs Basilisk within your project directory
- Sets up environment variables locally (in `.project_config`)
- Automatically detects your OS (MacOS or Linux) and uses appropriate configuration
- Verifies the installation

#### Basic usage:
```bash
# Run the script to install or use existing installation
./reset_install_requirements.sh

# For a fresh installation (removes existing one if present)
./reset_install_requirements.sh --hard

# Load the environment settings for your current shell session
source .project_config
```

The script will create a `.project_config` file in your project root with the necessary environment variables. This approach avoids modifying your global shell configuration and keeps the Basilisk installation contained within your project.

### Windows Subsystem for Linux (WSL) Compatibility

Testing on WSL is currently incomplete. In principle, the Linux installation instructions should work for WSL environments. If you encounter any issues while installing or running Basilisk on WSL, please report them by [opening a bug report](https://github.com/comphy-lab/Basilisk-101/issues/new?template=bug_report.md&labels=bug,wsl).

### Complete Installation Instructions

For more detailed installation instructions, including configuration for different systems, setting up environment variables, installing additional dependencies, and optional libraries, please refer to the official installation guide at [http://basilisk.fr/src/INSTALL](http://basilisk.fr/src/INSTALL).

### Running the codes

To use the make file do:
```bash
CFLAGS=-DDISPLAY=-1 make NAME-of-File.tst
```

## Reporting Issues and Feature Requests

We use GitHub Issues to track bugs, feature requests, and example requests for this course. When creating an issue, please select the appropriate template to help us address your needs efficiently.

### Issue Templates

#### Bug Report:
[Report here](https://github.com/comphy-lab/Basilisk-101/issues/new?template=bug_report.md)

- For problems with installation, compilation, or running code. 
Please include:
- Detailed description of the issue
- Your environment (OS, compiler version)
- Steps to reproduce
- Expected vs. actual behavior
- Error messages and logs
- Code snippets or files that demonstrate the issue

#### Feature/Topic Request:
[Report here](https://github.com/comphy-lab/Basilisk-101/issues/new?template=feature_request.md)
- For requesting specific topics or functionality
- Coverage of specific topics in the course
- New examples or tutorials
- Additional functionality in the codebase
- Improvements to existing materials

#### Example Request:
[Report here](https://github.com/comphy-lab/Basilisk-101/issues/new?template=example_request.md)
- For requesting specific examples that demonstrate:
- Particular Basilisk features
- Solutions to common problems
- Implementation of specific physics or numerical methods

#### General Question:
[Report here](https://github.com/comphy-lab/Basilisk-101/issues/new?template=general_question.md)
- For any other questions

### How to Create an Issue

1. Go to the ["Issues" tab](https://github.com/comphy-lab/Basilisk-101/issues) in the GitHub repository
2. Click the ["New Issue"](https://github.com/comphy-lab/Basilisk-101/issues/new/choose) button
3. Select the appropriate template from the options
4. Fill in the required information according to the template
5. Add relevant labels if available
6. Submit the issue