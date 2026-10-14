# Documentation style

Follow the [PyEnsembl writing guide](https://github.com/openvax/pyensembl/blob/edd183d/docs/dev/documentation-style.md)
and the [mhctools example](https://github.com/openvax/mhctools/tree/master/docs).

Put installation and the first useful example on the home page. Give complete
imports and concrete inputs, explain the result, and state prerequisites before
the reader runs a command. Show output only when it has been checked. Label
replaceable paths and templates.

Order guides by common tasks before advanced setup and implementation details.
Keep desktop sidebar groups expanded initially and individually collapsible.
Split a page when a substantial task or reference needs its own explanation.
Preserve useful old links and anchors when moving content.

Use ordinary text for concepts and library names; reserve inline code for
literal commands, API names, parameters, paths and values. Keep qualifications,
experimental status and scientific limitations beside the claims they qualify.

Use system fonts, a readable line length and restrained headings. Let wide
tables and code blocks scroll inside their containers. Avoid repetitive
navigation lists in prose and excessive callouts.

Install requirements-docs.txt and run ./docs.sh before opening a PR. Review
the home page and a reference page at desktop and narrow widths. Run the
repository's required lint and test commands and execute complete examples.
Documentation dependencies are separate from runtime dependencies.
