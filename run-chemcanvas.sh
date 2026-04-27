#!/bin/bash

# ChemCanvas Launcher Script for Fedora
# This script compiles UI resources and launches ChemCanvas without installation

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo -e "${GREEN}=== ChemCanvas Launcher ===${NC}\n"

# Check if Python 3 is installed
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}Error: Python 3 is not installed${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Python 3 found${NC}"

# Check for required dependencies
MISSING_DEPS=0

if ! python3 -c "import PyQt5" 2>/dev/null; then
    echo -e "${YELLOW}Missing: python3-pyqt5${NC}"
    MISSING_DEPS=1
fi

if ! command -v pyrcc5 &> /dev/null; then
    echo -e "${YELLOW}Missing: pyqt5-dev-tools${NC}"
    MISSING_DEPS=1
fi

if ! command -v pyuic5 &> /dev/null; then
    echo -e "${YELLOW}Missing: pyqt5-dev-tools${NC}"
    MISSING_DEPS=1
fi

#if [ $MISSING_DEPS -eq 1 ]; then
#    echo -e "${YELLOW}\nInstalling missing dependencies...${NC}"
#    echo "Running: sudo dnf install python3-pyqt5 pyqt5-devel"
#    sudo dnf install -y python3-pyqt5 pyqt5-devel
#    echo -e "${GREEN}Dependencies installed${NC}\n"
#fi

# Compile resources
echo -e "${YELLOW}Compiling resources...${NC}"

if [ ! -f "data/resources.qrc" ]; then
    echo -e "${RED}Error: data/resources.qrc not found${NC}"
    exit 1
fi

if [ ! -f "data/mainwindow.ui" ]; then
    echo -e "${RED}Error: data/mainwindow.ui not found${NC}"
    exit 1
fi

# Compile resource file
pyrcc5 -o chemcanvas/resources_rc.py data/resources.qrc
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Resources compiled${NC}"
else
    echo -e "${RED}Error compiling resources${NC}"
    exit 1
fi

# Compile UI file
pyuic5 -o chemcanvas/ui_mainwindow.py data/mainwindow.ui
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ UI compiled${NC}"
else
    echo -e "${RED}Error compiling UI${NC}"
    exit 1
fi

# Launch ChemCanvas
echo -e "${GREEN}\n✓ All resources compiled successfully${NC}"
echo -e "${GREEN}Launching ChemCanvas...${NC}\n"

python3 -m chemcanvas.main "$@"
