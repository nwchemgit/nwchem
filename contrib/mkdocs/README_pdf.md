
MKDOCS_SERVE=B ./build_website.sh

cd site

wkhtmltopdf --javascript-delay 45000 --enable-local-file-access print_page.html print_pdf.pdf
