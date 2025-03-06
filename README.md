# Basilisk-101

## Installation

Basilisk requires a C99-compliant compiler and GNU make. Installation can be done in two ways:

## Global installation

### Using darcs (recommended)
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

### Using a tarball
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
### Add to your shell configuration (.bashrc or .zshrc)
```bash
echo "export BASILISK=$PWD" >> ~/.bashrc
echo 'export PATH=$PATH:$BASILISK' >> ~/.bashrc
```
Or for zsh users:
```bash
echo "export BASILISK=$PWD" >> ~/.zshrc
echo 'export PATH=$PATH:$BASILISK' >> ~/.zshrc
```

## Repository level installation

For project-specific installations, you can use the provided `reset_install_requirements.sh` script which:
- Installs Basilisk within your project directory
- Sets up environment variables locally (in `.project_config`)
- Automatically detects your OS (MacOS or Linux) and uses appropriate configuration (config.osx for Mac, config.gcc for Linux)
- Verifies the installation

### Basic usage:
```bash
# Run the script to install or use existing installation
./reset_install_requirements.sh

# For a fresh installation (removes existing one if present)
./reset_install_requirements.sh --hard

# Load the environment settings for your current shell session
source .project_config
```

The script will create a `.project_config` file in your project root with the necessary environment variables. This approach avoids modifying your global shell configuration and keeps the Basilisk installation contained within your project.

## Windows Subsystem for Linux (WSL) Compatibility

Testing on WSL is currently incomplete. In principle, the Linux installation instructions should work for WSL environments. If you encounter any issues while installing or running Basilisk on WSL, please report them by opening an issue on GitHub.

## Complete installation instructions
For complete installation instructions, including:
- Configuration for different systems
- Setting up environment variables
- Installing additional dependencies
- Optional libraries

Please refer to the official installation guide at:
[http://basilisk.fr/src/INSTALL](http://basilisk.fr/src/INSTALL)