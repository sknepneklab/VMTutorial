from skbuild import setup

setup(
    name="vmtutorial",
    version="0.5",
    description="Vertex Model Tutorial",
    author="Rastko Sknepnek",
    license="MIT",
    packages=["vmtutorial"],
    package_dir={"": "src"},
    cmake_install_dir="src/vm",
    python_requires=">=3.8",
)
