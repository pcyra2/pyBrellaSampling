from setuptools import setup
# set up using "pip install -e ."
setup(
    name='pyBrellaSampling',
    version='2.0',
    py_modules=['pyBrellaSampling'],
    entry_points={
        'console_scripts': [
            'pyBrella = pyBrellaSampling.pyBrella:main',
            'Standalone = pyBrellaSampling.Standalone:main',
            'Benchmark = pyBrellaSampling.Benchmark:main',
            'convert = pyBrellaSampling.Tools.utils:convert',
            'pyscfQMMM = pyBrellaSampling.Code.Tools.QM.run_pyscf:main',
            'qmmm_sniffer = pyBrellaSampling.Code.Tools.QM.qmmm_sniffer:main',
            'dm21_test = pyBrellaSampling.Code.Tools.QM.testDM21:main'
        ],
    },
)
