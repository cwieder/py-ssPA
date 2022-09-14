
::: sspa.sspa_ora
    handler: python
    options:
      show_root_heading: false
      show_source: false

::: sspa.identifier_conversion.identifier_conversion
    handler: python
    options:
      members:
        - identifier_conversion
        - map_identifiers
      show_root_heading: false
      show_source: false

::: sspa.identifier_conversion.map_identifiers
    handler: python
    options:
      members:
        - map_identifiers
      show_root_heading: false
      show_source: false

::: sspa.download_pathways
    handler: python
    options:
      members:
        - download_reactome
        - download_KEGG
        - MetExplorePaths
      show_root_heading: true
      show_source: false