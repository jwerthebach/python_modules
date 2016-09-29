from setuptools import setup, find_packages

setup(
	name = "jopymods",
	version = "0.1",
	description = "Module with python functions for analysis",
	url = "https://github.com/jwerthebach/python_modules",
	author = "Johannes Werthebach",
	author_email = "johannes.werthebach@tu-dortmund.de",
	license = "MIT",
	install_requires = [
		"numpy",
		"scipy",
		"matplotlib",
		"pandas",
		"tables",
		"tqdm",
		"seaborn",
		],
	packages = find_packages(),
)