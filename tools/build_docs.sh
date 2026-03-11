#!/bin/bash
#
# This file is part of `pyfraglib`, a software suite to calculate fragmentomics
# features from cfDNA and perform downstream analyses.
#
# Copyright (C) 2025 Daniel Schütte, daniel.schuette@iccb-cologne.org
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version. This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details. You should have received a copy of the GNU General Public
# License along with this program. If not, see <https://www.gnu.org/licenses/>.
RED="\033[0;31m"
GREEN="\033[0;32m"
YELLOW="\033[1;33m"
NC="\033[0m" # no Color

set -e

setup() {
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    DOCS_DIR="$(cd "$SCRIPT_DIR/../docs" && pwd)"
    cd "$DOCS_DIR"

    echo -e "${GREEN}pyfraglib Documentation Builder${NC}"
    if [ ! -f "source/conf.py" ]; then
        echo -e "${RED}Error: conf.py not found. Are you in the docs directory?${NC}"
        exit 1
    fi

    echo -e "${YELLOW}Checking dependencies...${NC}"
    if ! python -c "import sphinx" 2>/dev/null; then
        echo -e "${YELLOW}Installing documentation dependencies...${NC}"
        pip install -r requirements.txt
    fi
}

build_html() {
    echo -e "${YELLOW}Building HTML documentation...${NC}"
    mkdir -p build

    echo -e "${YELLOW}Running sphinx-build...${NC}"
    sphinx-build -b html source build/html
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}HTML documentation built successfully.${NC}"
        echo -e "${GREEN}HTML files: docs/build/html/index.html${NC}"
        echo -e "${GREEN}Open docs/build/html/index.html in your browser.${NC}"
    else
        echo -e "${RED}HTML documentation build failed.${NC}"
        exit 1
    fi
}

build_pdf() {
    echo -e "${YELLOW}Building PDF documentation...${NC}"
    mkdir -p build

    if ! command -v pdflatex &> /dev/null; then
        echo -e "${RED}Error: pdflatex not found. Please install LaTeX.${NC}"
        exit 1
    fi

    echo -e "${YELLOW}Running sphinx-build for LaTeX...${NC}"
    sphinx-build -b latex source build/latex
    if [ $? -ne 0 ]; then
        echo -e "${RED}LaTeX build failed.${NC}"
        exit 1
    fi

    echo -e "${YELLOW}Running pdflatex...${NC}"
    cd build/latex
    set +e
    for i in {1..3}; do
        echo -e "${YELLOW}pdflatex run $i/3...${NC}"
        pdflatex -interaction=nonstopmode pyfraglib_documentation.tex
        if [ ! -f "pyfraglib_documentation.pdf" ]; then
            echo -e "${RED}pdflatex failed on run $i - no PDF generated.${NC}"
            exit 1
        fi
        if [ $i -eq 1 ] && command -v makeindex &> /dev/null; then
            echo -e "${YELLOW}Running makeindex...${NC}"
            makeindex -s python.ist pyfraglib_documentation.idx 2>/dev/null || true
        fi
    done
    set -e

    cd "$DOCS_DIR"
    if [ -f "build/latex/pyfraglib_documentation.pdf" ]; then
        cp build/latex/pyfraglib_documentation.pdf pyfraglib_documentation.pdf
        echo -e "${GREEN}PDF documentation built successfully.${NC}"
        echo -e "${GREEN}PDF file: docs/pyfraglib_documentation.pdf${NC}"
    else
        echo -e "${RED}PDF generation failed.${NC}"
        exit 1
    fi
}

case "${1}" in
    html)
        setup
        build_html
        ;;
    pdf)
        setup
        build_pdf
        ;;
    *)
        echo "Usage: $0 [html|pdf]"
        echo "  html  - Build HTML documentation (default)"
        echo "  pdf   - Build PDF documentation"
        exit 1
        ;;
esac
