import re
line = "     1                 -87.6136               0.00411523"
line = line.strip()
_match = re.match(r"^([0-9]+)(\s+)(-?[0-9]+\.[0-9]+)(\s+)([01](\.[0-9]+)?)(.*)", line)
print(_match.groups())
line = "     1                  -91.611                        1                 -91.4246                        1"
line = line.strip()
_match = re.match(r"^([0-9]+)(\s+)(-?[0-9]+\.[0-9]+)(\s+)([01](\.[0-9]+)?)(.*)", line)
print(_match.groups())
_match = re.match(r"^(-?[0-9]+\.[0-9]+)(\s+)([01](\.[0-9]+)?)(.*)", _match.group(7).strip())
print(_match.groups())