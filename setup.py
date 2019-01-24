from setuptools import setup

setup(name='claws',
      version='0.1',
      description='Framework to handle the Claws data',
#      url='http://github.com/storborg/funniest',
      author='Hendrik Windel',
      author_email='hwindel@mpp.mpg.de',
      license='MIT',
      packages=['claws'],
      install_requires=['pandas', 'configparser'],
      zip_safe=False)

