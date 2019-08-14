# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)


def stringify_keys(d):
    # taken from https://stackoverflow.com/a/51051641
    """Convert a dict's keys to strings if they are not."""
    keys = list(d.keys())
    for key in keys:

        # check inner dict
        if isinstance(d[key], dict):
            value = stringify_keys(d[key])
        else:
            value = d[key]

        # convert nonstring to string if needed
        if not isinstance(key, str):
            try:
                d[str(key)] = value
            except Exception:
                try:
                    d[repr(key)] = value
                except Exception:
                    raise

            # delete old key
            d.pop(key, None)
    return d
