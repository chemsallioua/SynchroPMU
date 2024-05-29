from setuptools import setup, find_packages

setup(
    name="pmu-estimator",
    version="1.2.0",
    description="Python bindings for PmuEstimator",
    author="Chemseddine Allioua",
    author_email="chemesallioua@gmail.com",
    url="https://github.com/chemsallioua/PmuEstimator",
    packages=find_packages(),
    py_modules=["pmu_estimator"],
    classifiers=[
        "Development Status :: 1 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires='>=3.6',
)
