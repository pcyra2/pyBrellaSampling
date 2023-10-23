from setuptools import setup
# set up using "pip install -e ."
setup(
    name='pyBrellaSampling',
    version='1.0',
    py_modules=['pyBrellaSampling'],
    entry_points={
        'console_scripts': [
            'pyBrella = pyBrellaSampling.pyBrella:main',
            'Standalone = pyBrellaSampling.Standalone:main_cli',
        ],
    },
)
