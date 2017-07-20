from skbuild import setup

setup(
        name='goofit',
        version='2.1',
        description='GooFit fitting package',
        author='Henry Schreiner',
        author_email='hschrein@cern.ch',
        url='https://goofit.github.io',
        platforms = ["POSIX"],
        provides = ["goofit"],
        license="LGPL 3.0",
        packages=['goofit']
        # description="A parallel fitting package."
        )

