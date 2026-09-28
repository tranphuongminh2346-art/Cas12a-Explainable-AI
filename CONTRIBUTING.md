# Contributing to cas12a-xai

We welcome contributions from computational biologists, bioinformaticians, and machine learning researchers!

## How to Contribute

### 1. Reporting Bugs & Asking Questions
- Search the [Issues](https://github.com/tranphuongminh2346-art/Cas12a-Explainable-AI/issues) tracker to see if the issue has already been reported.
- If not, open a new issue describing:
  - The operating system and Python version
  - Exact command or code snippet executed
  - Full traceback and error output

### 2. Submitting Pull Requests
1. Fork the repository on GitHub.
2. Clone your fork locally:
   ```bash
   git clone https://github.com/<your-username>/Cas12a-Explainable-AI.git
   cd Cas12a-Explainable-AI
   ```
3. Create a feature branch:
   ```bash
   git checkout -b feature/new-biophysical-descriptor
   ```
4. Install editable development mode with testing tools:
   ```bash
   pip install -e .
   pip install pytest flake8
   ```
5. Run the test suite before making changes:
   ```bash
   python -m unittest discover tests
   ```
6. Implement your changes, write tests, and ensure code styling adheres to PEP 8.
7. Commit and push your branch:
   ```bash
   git commit -m "feat: add RNA secondary structure folding stability"
   git push origin feature/new-biophysical-descriptor
   ```
8. Open a Pull Request with a clear description of the biological rationale and benchmarking results.

## Code of Conduct
Please be respectful and constructive in all code reviews, discussions, and issues.
