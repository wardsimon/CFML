#!/bin/bash
set -e -u -x

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" --plat "$PLAT" -w /io/dist/
    fi
}


# Install a system package required by our library
#yum install -y atlas-devel
# yum install -y libaec-dev libz-dev libsz2


# Compile wheels
for PYBIN in /opt/python/*/bin; do
#    "${PYBIN}/pip" install -r /io/dev-requirements.txt
    "${PYBIN}/pip" wheel /io/ --no-deps --use-feature=in-tree-build -w /io/dist/
    rm -rf /io/build/
done

# Bundle external shared libraries into the wheels
for whl in /io/dist/*.whl; do
    repair_wheel "$whl"
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip" install CFML --no-index -f /io/dist
#    (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
done