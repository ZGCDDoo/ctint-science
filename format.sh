#!/bin/bash

# inspired by  https://github.com/gabime/spdlog/blob/v1.x/format.sh

echo -n "Running dos2unix     "
find . -name "*\.hpp" -o -name "*\.cpp"|grep -v bundled|xargs -I {} sh -c "dos2unix '{}' 2>/dev/null; echo -n '.'"
echo
echo -n "Running clang-format "
find . -name "*\.hpp" -not -path "./deps/*" -o -name "*\.cpp"|grep -v bundled|xargs -I {} sh -c "clang-format -i {}; echo -n '.'"
echo
