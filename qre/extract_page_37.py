"""Extract page 37 and surrounding pages which contain equation 145."""

import PyPDF2

pdf_path = "/vast/home/bkkrueger/.claude/projects/-vast-home-bkkrueger-qc4ws-lisdi-git-qhat-worktrees-bkk-commutator-qre/1b372140-31cc-4899-b3c0-378f66f8c5be/tool-results/webfetch-1776705006385-jfcs70.pdf"

with open(pdf_path, 'rb') as file:
    reader = PyPDF2.PdfReader(file)

    # Section 5.1 is listed as page 37 in the table of contents
    # PDF pages are 0-indexed, but the paper page numbers start later
    # Let's search for pages containing "5.1" or "First- and second-order"

    print("Searching for Section 5.1 (First- and second-order error bounds)...\n")

    for page_num in range(35, 45):  # Check pages around 37
        page = reader.pages[page_num]
        text = page.extract_text()

        if '5.1' in text or 'First- and second-order' in text or 'Equation 145' in text or '(145)' in text:
            print("="*70)
            print(f"PAGE {page_num + 1} (0-indexed: {page_num})")
            print("="*70)
            print(text)
            print("\n\n")
