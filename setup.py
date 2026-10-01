from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="cas12a-xai",
    version="0.1.2",
    author="Minh Tran",
    author_email="tranphuongminh2346@gmail.com",
    description="Explainable Machine Learning for CRISPR-Cas12a On-Target Cleavage Efficiency",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tranphuongminh2346-art/Cas12a-Explainable-AI",
    project_urls={
        "Bug Tracker": "https://github.com/tranphuongminh2346-art/Cas12a-Explainable-AI/issues",
        "Documentation": "https://github.com/tranphuongminh2346-art/Cas12a-Explainable-AI#readme",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    package_data={"cas12a_xai": ["models/*.txt"]},
    include_package_data=True,
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.22.0",
        "pandas>=1.4.0",
        "scipy>=1.8.0",
        "scikit-learn>=1.0.0",
        "lightgbm>=3.3.0",
        "shap>=0.41.0",
    ],
    entry_points={
        "console_scripts": [
            "cas12a-xai=cas12a_xai.cli:main",
        ],
    },
)
