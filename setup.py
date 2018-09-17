from setuptools import setup

packages = [
				'gb_code',
       ]

INSTALL_REQUIRES = (
    ['numpy >= 1.14.0'] )

setup(
	name='GB_code',
    python_requires='>3.5.1',
    version='1.0.0',
    author='R.Hadian',
    author_email='shahrzadhadian@gmail.com',
    packages=packages,
				install_requires=INSTALL_REQUIRES,
				entry_points = {
        'console_scripts': [
            'csl_generator = gb_code.csl_generator:main',
            'gb_generator = gb_code.gb_generator:main',
        ],
    },
			 url='https://github.com/oekosheri/GB_code',
    license='LICENSE',
    description='A GB generation code',
    long_description=open('README.md').read(),
)
