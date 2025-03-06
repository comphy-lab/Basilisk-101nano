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

Testing on WSL is currently incomplete. In principle, the Linux installation instructions should work for WSL environments. If you encounter any issues while installing or running Basilisk on WSL, please report them by [opening a bug report](https://github.com/USERNAME/REPOSITORY/issues/new?template=bug_report.md&labels=bug,wsl).

## Reporting Issues and Feature Requests

We use GitHub Issues to track bugs, feature requests, and example requests for this course. When creating an issue, please select the appropriate template to help us address your needs efficiently.

### Issue Templates

We provide the following templates for different types of issues:

#### [Bug Report Template](https://github.com/USERNAME/REPOSITORY/issues/new?template=bug_report.md)
Use this template when you encounter problems with installation, compilation, or running Basilisk code. Include:
- Detailed description of the issue
- Your environment (OS, compiler version)
- Steps to reproduce
- Expected vs. actual behavior
- Error messages and logs
- Code snippets or files that demonstrate the issue

#### [Feature/Topic Request Template](https://github.com/USERNAME/REPOSITORY/issues/new?template=feature_request.md)
Use this when you'd like to request:
- Coverage of specific topics in the course
- New examples or tutorials
- Additional functionality in the codebase
- Improvements to existing materials

#### [Example Request Template](https://github.com/USERNAME/REPOSITORY/issues/new?template=example_request.md)
Use this template when requesting specific examples that demonstrate:
- Particular Basilisk features
- Solutions to common problems
- Implementation of specific physics or numerical methods

#### [General Question Template](https://github.com/USERNAME/REPOSITORY/issues/new?template=general_question.md)
For any questions about Basilisk usage, the course material, or computational fluid dynamics concepts that aren't covered by the other templates.

### How to Create an Issue

1. Go to the ["Issues" tab](https://github.com/USERNAME/REPOSITORY/issues) in the GitHub repository
2. Click the ["New Issue"](https://github.com/USERNAME/REPOSITORY/issues/new/choose) button
3. Select the appropriate template from the options
4. Fill in the required information according to the template
5. Add relevant labels if available
6. Submit the issue

You can also create issues directly using these links:
- [Report a bug](https://github.com/USERNAME/REPOSITORY/issues/new?template=bug_report.md)
- [Request a feature or topic](https://github.com/USERNAME/REPOSITORY/issues/new?template=feature_request.md)
- [Request an example](https://github.com/USERNAME/REPOSITORY/issues/new?template=example_request.md)
- [Ask a question](https://github.com/USERNAME/REPOSITORY/issues/new?template=general_question.md)
- [Submit another type of issue](https://github.com/USERNAME/REPOSITORY/issues/new?template=blank_issue.md)

Please provide as much relevant information as possible to help us understand and address your issue efficiently.

## Complete installation instructions
For complete installation instructions, including:
- Configuration for different systems
- Setting up environment variables
- Installing additional dependencies
- Optional libraries

Please refer to the official installation guide at:
[http://basilisk.fr/src/INSTALL](http://basilisk.fr/src/INSTALL)