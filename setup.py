# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import re

from setuptools import setup, find_packages



package_name = "gtfparse"
current_directory = os.path.dirname(__file__)
readme_filename = 'README.md'
readme_path = os.path.join(current_directory, readme_filename)
github_url = "https://github.com/openvax/%s" % package_name


readme_markdown = ""
try:
    with open(readme_path, 'r') as f:
        readme_markdown = f.read()
except Exception as e:
    print(e)
    print("Failed to open %s" % readme_path)

with open('gtfparse/version.py', 'r') as f:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        f.read(),
        re.MULTILINE).group(1)

with open("requirements.txt") as f:
    requirements = [l.strip() for l in f]


if __name__ == '__main__':
    setup(
        name=package_name,
        packages=find_packages(),
        version=version,
        description="GTF Parsing",
        long_description=readme_markdown,
        long_description_content_type='text/markdown',
        url=github_url,
        author="Alex Rubinsteyn",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires=requirements,
        package_data={
            package_name: ['../requirements.txt'],
        },
    )
