#! /bin/bash

echo "testing container..."

# make sure current dir is repo dir
cd $( dirname "${BASH_SOURCE[0]}" )

TEST_IMAGE="facets_ch_test_image"

if [ "$1" = "--skip-build" ]; then
    echo "skipping build..."
else
    echo "building image, to skip run with --skip-build..."
    docker build -t $TEST_IMAGE .
fi

# see https://explainshell.com/explain?cmd=set+-euxo%20pipefail
set -euxo pipefail

# run tox inside the container
echo "testing docker image..."
docker run --rm $TEST_IMAGE --help
docker run --rm -it --entrypoint "" -v `pwd`:/test -w /test $TEST_IMAGE bash /test/tests/run_test.sh
echo "tests finished..."
