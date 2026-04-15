import os
import re
import sys
import sysconfig
import platform
import subprocess

from pathlib import Path
from packaging.version import Version
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

ROOT_DIR = Path(__file__).resolve().parent
README_TEXT = (ROOT_DIR / "README.md").read_text(encoding="utf-8")

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extension: "
                               + ", ".join(e.name for e in self.extensions))

        cmake_version = Version(re.search(r"version\s*([\d.]+)", out.decode()).group(1))
        if platform.system() == "Windows" and cmake_version < Version("3.14"):
            raise RuntimeError("CMake >= 3.14 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        prefix = sysconfig.get_config_var("LIBDIR")

        # CMake arguments
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
            "-DPYTHON_LIBRARY_DIR={}".format(prefix)
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        # Platform-specific settings
        if platform.system() == "Windows":
            cmake_args += ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            # Linux/Unix settings
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j2"]

        # Add environment variables for compilation
        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get("CXXFLAGS", ""), self.distribution.get_version()).replace('"', '\\"')

        # Create the build directory if it does not exist
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # Run CMake
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)

        # Move the generated output to the correct location
        self.move_output(ext)

    def move_output(self, ext):
        extdir = Path(self.build_lib).resolve()
        dest_path = Path(self.get_ext_fullpath(ext.name)).resolve()
        source_path = extdir / self.get_ext_filename(ext.name)
        dest_directory = dest_path.parent

        # Check whether the file was generated
        if not source_path.exists():
            raise RuntimeError(
                f"Output file {source_path} not found. Check whether CMake generated the .so file correctly."
            )

        # Create the destination directory if needed and move the file
        dest_directory.mkdir(parents=True, exist_ok=True)
        self.copy_file(source_path, dest_path)

        



        
system = platform.system()

if system == "Linux":
    native_ext = "*.so"
elif system == "Darwin":
    native_ext = "*.so"
elif system == "Windows":
    native_ext = "*.dll"
else:
    raise RuntimeError(f"Unsupported platform: {system}!")

setup(
    name="morphoTreeAdjust",
    version="0.2",
    description="Core C++/Python implementation of MorphoTreeAdjust for dynamic component-tree adjustment.",
    long_description=README_TEXT,
    long_description_content_type="text/markdown",
    author="Wonder Alexandre Luz Alves",
    author_email="worderalexandre@gmail.com",
    license="GPL-3.0",
    url="https://github.com/wonderalexandre/MorphoTreeAdjust",
    keywords="machine learning, morphological trees, mathematical morphology",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Image Processing",
        "Programming Language :: Python",
        "Programming Language :: C++",
    ],
    packages=find_packages(where="python", include=["morphoTreeAdjust", "morphoTreeAdjust.*"]),
    package_dir={"": "python"},
    ext_modules=[CMakeExtension('morphoTreeAdjust.morphoTreeAdjust')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    include_package_data=True,
    package_data={"morphoTreeAdjust": ["*.py", native_ext]},
    extras_require={
        "viz": ["numpy", "matplotlib"],
    }
)
