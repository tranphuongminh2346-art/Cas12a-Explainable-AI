"""
One-click standalone script to patch all Jupyter notebooks in notebooks/
Ensures axis labels match 'Predicted Guide Efficiency (Arbitrary / Indel Score)',
harmonizes author affiliation to Independent Researcher via regex,
and annotates poly_t_terminator as 23-bp target window (n=1,042; in-guide n=927).
"""
import os
import re

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NB_DIR = os.path.join(BASE_DIR, "notebooks")

def patch_all_notebooks():
    print(f"[*] Scanning notebooks in {NB_DIR}...")
    if not os.path.exists(NB_DIR):
        print(f"[!] Directory not found: {NB_DIR}")
        return
    for fname in os.listdir(NB_DIR):
        if not fname.endswith(".ipynb"):
            continue
        fpath = os.path.join(NB_DIR, fname)
        try:
            with open(fpath, "r", encoding="utf-8") as f:
                content = f.read()
            original_content = content
            
            # Harmonize axis labels
            content = content.replace(
                'axes[0, 0].set_xlabel("Predicted On-Target Efficiency (%)")',
                'axes[0, 0].set_xlabel("Predicted Guide Efficiency (Arbitrary / Indel Score)")'
            )
            content = content.replace(
                'axes[0, 0].set_xlabel(\\"Predicted On-Target Efficiency (%)\\")',
                'axes[0, 0].set_xlabel(\\"Predicted Guide Efficiency (Arbitrary / Indel Score)\\")'
            )
            
            # Regex-based author affiliation standardization
            content = re.sub(
                r'(\*\*Author\*\*:?\s*Minh Tran\s*\()[^)]*(\))',
                r'\1Independent Researcher, Perth, Western Australia, Australia\2',
                content
            )
            content = re.sub(
                r'(Author:\s*Minh Tran\s*\()[^)]*(\))',
                r'\1Independent Researcher, Perth, Western Australia, Australia\2',
                content
            )
            
            # Harmonize target journal
            content = re.sub(
                r'(\*\*Target Journal\*\*:?\s*)[^\n\\]+',
                r'\1Computational Biology & Bioinformatics Journal (Traditional / Subscription Route)',
                content
            )
            
            if "'poly_t_terminator': int('TTTT' in sp)" in content and "# 23-bp target window" not in content:
                content = content.replace(
                    "'poly_t_terminator': int('TTTT' in sp)",
                    "'poly_t_terminator': int('TTTT' in sp)  # 23-bp target window (n=1,042; in-guide n=927 drives -38.38% deficit)"
                )
                
            if content != original_content:
                with open(fpath, "w", encoding="utf-8") as f:
                    f.write(content)
                print(f"  [OK] Successfully patched: {fname}")
            else:
                print(f"  [-] Already clean: {fname}")
        except Exception as e:
            print(f"  [!] Error patching {fname}: {e}")

if __name__ == "__main__":
    patch_all_notebooks()
