"""
One-click standalone script to patch all Jupyter notebooks in notebooks/
Ensures axis labels match 'Predicted Guide Efficiency (Arbitrary / Indel Score)'
and annotates poly_t_terminator as 23-bp target window (n=1,042; in-guide n=927).
"""
import os

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
            modified = False
            axis_replacements = [
                (
                    'axes[0, 0].set_xlabel("Predicted On-Target Efficiency (%)")',
                    'axes[0, 0].set_xlabel("Predicted Guide Efficiency (Arbitrary / Indel Score)")',
                ),
                (
                    'axes[0, 0].set_xlabel(\\"Predicted On-Target Efficiency (%)\\")',
                    'axes[0, 0].set_xlabel(\\"Predicted Guide Efficiency (Arbitrary / Indel Score)\\")',
                ),
                (
                    'School of Molecular Sciences & Centre for Applied Bioinformatics, UWA / Independent Researcher',
                    'Independent Researcher, Perth, Western Australia, Australia',
                ),
                (
                    'University of Western Australia / Independent Researcher',
                    'Independent Researcher, Perth, Western Australia, Australia',
                ),
                (
                    'School of Molecular Sciences, The University of Western Australia',
                    'Independent Researcher, Perth, Western Australia, Australia',
                ),
                (
                    'School of Molecular Sciences, University of Western Australia',
                    'Independent Researcher, Perth, Western Australia, Australia',
                ),
                (
                    'Minh Tran (School of Molecular Sciences, University of Western Australia)',
                    'Minh Tran (Independent Researcher, Perth, Western Australia, Australia)',
                ),
                (
                    'Minh Tran (Independent Researcher)',
                    'Minh Tran (Independent Researcher, Perth, Western Australia, Australia)',
                ),
            ]
            for old, new in axis_replacements:
                if old in content:
                    content = content.replace(old, new)
                    modified = True
            if "'poly_t_terminator': int('TTTT' in sp)" in content and "# 23-bp target window" not in content:
                content = content.replace(
                    "'poly_t_terminator': int('TTTT' in sp)",
                    "'poly_t_terminator': int('TTTT' in sp)  # 23-bp target window (n=1,042; in-guide n=927 drives -38.38% deficit)"
                )
                modified = True
            if modified:
                with open(fpath, "w", encoding="utf-8") as f:
                    f.write(content)
                print(f"  [OK] Successfully patched: {fname}")
            else:
                print(f"  [-] Already clean: {fname}")
        except Exception as e:
            print(f"  [!] Error patching {fname}: {e}")

if __name__ == "__main__":
    patch_all_notebooks()
