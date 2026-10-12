# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""YAML a user gives foldmap (layout files, names files, style files), read safely: plain values only, no
anchors (&a) or aliases (*a). An alias repeats a node without copying text, so a few hundred bytes of nested
aliases ("billion laughs") expand to gigabytes once read; none of foldmap's files needs one."""

from __future__ import annotations

import yaml


class _Refused(Exception):
    """An anchor or alias in the document (raised from inside the YAML reader)."""


class _NoAliases(yaml.SafeLoader):
    """SafeLoader that stops at the first anchor or alias, before any node is built from it."""

    def compose_node(self, parent, index):
        event = self.peek_event()
        if isinstance(event, yaml.AliasEvent) or getattr(event, "anchor", None) is not None:
            raise _Refused
        return super().compose_node(parent, index)


def read_yaml(text: str, kind: str, where: str):
    """The YAML document in `text` (None when empty), read with no anchors or aliases. kind names the file in the
    messages ("layout", "names", "style"), where names the file itself. Every failure is a ValueError saying what
    is wrong: an anchor or alias, invalid YAML, or brackets nested too deeply for the reader."""
    loader = _NoAliases(text)
    try:
        return loader.get_single_data()
    except _Refused:
        raise ValueError(f"{where}: aliases are not allowed in a {kind} file (no &anchor or *alias)") from None
    except yaml.YAMLError as err:
        raise ValueError(f"{where}: not valid YAML: {err}") from None
    except RecursionError:  # brackets nested thousands deep: the YAML reader recurses once per level
        raise ValueError(f"{where}: nested too deeply") from None
    except ValueError:  # the constructor's own limit, e.g. an integer of more than 4,300 digits
        raise ValueError(f"{where}: a value is too long to read (e.g. a number of thousands of digits)") from None
    finally:
        loader.dispose()
