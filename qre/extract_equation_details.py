"""Extract pages 38-40 to see the full second-order formula."""

import PyPDF2

pdf_path = "/vast/home/bkkrueger/.claude/projects/-vast-home-bkkrueger-qc4ws-lisdi-git-qhat-worktrees-bkk-commutator-qre/1b372140-31cc-4899-b3c0-378f66f8c5be/tool-results/webfetch-1776705006385-jfcs70.pdf"

with open(pdf_path, 'rb') as file:
    reader = PyPDF2.PdfReader(file)

    # Pages 37-40 (0-indexed 36-39) should have the full analysis
    for page_num in range(37, 40):
        page = reader.pages[page_num]
        text = page.extract_text()

        print("="*70)
        print(f"PAGE {page_num + 1} (0-indexed: {page_num})")
        print("="*70)
        print(text)
        print("\n\n")
