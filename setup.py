import sys
import os
import shutil
import subprocess

from setuptools import setup
from setuptools.command.build_ext import build_ext

class BuildExt(build_ext):
    def run(self):
        cmd = ["make", "-j", "20"]
        if subprocess.call(cmd) != 0:
            sys.exit(-1)
        shutil.move("raspa3/raspa3.so", os.path.join(self.build_lib, "raspa3"))
        super().run()


setup(
    name="raspa3",
    version="3.0.0",
    packages=["raspa3"],
    has_ext_modules= lambda: True,
    cmdclass={"build_ext": BuildExt},
)