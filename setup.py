from skbuild import setup

setup(
    name="VMToolkit",
    version="0.1",
    description="Vertex Model Tutorial",
    author="Rastko Sknepnek",
    license="MIT",
    packages=['VMToolkit', 'VMToolkit.VM', 'VMToolkit.VMAnalysis', 'VMToolkit.VMAnalysis.Graner', 'VMToolkit.VMAnalysis.utils', 'VMToolkit.VMAnalysis.utils.HalfEdge'],
    cmake_install_dir='VMToolkit/VM',
    python_requires=">=3.8",
)