# Copyright (c) 2015. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import
import logging
from collections import OrderedDict

from six.moves import intern

from .util import memory_usage

def expand_attribute_strings(
        attribute_strings,
        quote_char='\"',
        missing_value=""):
    """
    The last column of GTF has a variable number of key value pairs
    of the format: "key1 value1; key2 value2;"
    Parse these into a dictionary mapping each key onto a list of values,
    where the value is None for any row where the key was missing.

    Parameters
    ----------
    attribute_strings : list of str

    quote_char : str
        Quote character to remove from values

    missing_value : any
        If an attribute is missing from a row, give it this value.

    Returns OrderedDict of column->value list mappings, in the order they
    appeared in the attribute strings.
    """
    logging.debug(
        "Memory usage before expanding GTF attributes: %0.4f MB" % (
            memory_usage(),))
    n = len(attribute_strings)

    extra_columns = {}
    column_order = []

    # Split the semi-colon separated attributes in the last column of a GTF
    # into a list of (key, value) pairs.
    kv_generator = (
        # We're slicing the first two elements out of split() because
        # Ensembl release 79 added values like:
        #   transcript_support_level "1 (assigned to previous version 5)";
        # ...which gets mangled by splitting on spaces.
        #
        # TODO: implement a proper parser!
        (i, kv.strip().split(" ", 2)[:2])
        for (i, attribute_string) in enumerate(attribute_strings)
        for kv in attribute_string.split(";")
        # need at least 3 chars for minimal entry like 'k v'
        if len(kv) > 2 and " " in kv
    )

    #
    # SOME NOTES ABOUT THE BIZARRE STRING INTERNING GOING ON BELOW
    #
    # While parsing millions of repeated strings (e.g. "gene_id" and "TP53"),
    # we can save a lot of memory by making sure there's only one string
    # object per unique string. The canonical way to do this is using
    # the 'intern' function. One problem is that Py2 won't let you intern
    # unicode objects, so to get around this we call intern(str(...)).
    #
    # It also turns out to be faster to check interned strings ourselves
    # using a local dictionary, hence the two dictionaries below
    # and pair of try/except blocks in the loop.
    column_interned_strings = {}
    value_interned_strings = {}

    for i, (column_name, value) in kv_generator:
        try:
            column_name = column_interned_strings[column_name]
            column = extra_columns[column_name]
        except KeyError:
            column_name = intern(str(column_name))
            column_interned_strings[column_name] = column_name
            column = [missing_value] * n
            extra_columns[column_name] = column
            column_order.append(column_name)

        value = value.replace(quote_char, "") if quote_char in value else value

        try:
            value = value_interned_strings[value]
        except KeyError:
            value = intern(str(value))
            value_interned_strings[value] = value

        column[i] = value

    logging.debug(
        "Memory usage after expanding GTF attributes: %0.4f MB" % (
            memory_usage(),))
    logging.info("Extracted GTF attributes: %s" % column_order)
    return OrderedDict(
        (column_name, extra_columns[column_name])
        for column_name in column_order)
