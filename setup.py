from setuptools import setup

packages = [ 
				'gb_code',
       ]


setup(
			 name='GB_code',
    version='0.1.0',
    author='R.Hadian',
    author_email='shahrzadhadian@gmail.com',
    packages=packages,
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
