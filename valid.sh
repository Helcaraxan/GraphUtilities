#!/bin/bash

#
# Script for building and validating the GraphUtilities library with different build options.
# -> Used for testing the library after modifications and before commiting
# -> Used during Continuous Integration to decide wether new commits produce a stable library
#

# Get the path to the place were the library resides and move to it
REPODIR="$(cd "$(dirname "${BASHSOURCE[0]}")" && pwd)"

cd "${REPODIR}"

# Remove any existing build directory
rm -rf "${REPODIR}/build"


# -------------------
# Sub-build procedure
# -------------------
runSubBuild() {
	# Create the sub-directory for this build and move there
	echo "Creating sub-dir for build: ${REPODIR}/build/$1"
	mkdir -p "${REPODIR}/build/$1"
	echo "Moving to : ${REPODIR}/build/$1"
	cd "${REPODIR}/build/$1"

	# Running cmake for this build with the specified options
	if [[ -z "$2" ]]; then
		echo "Running cmake with default options"
		cmake "${REPODIR}/." "-DCMAKE_BUILD_TYPE=Release"
	else
		echo "Running cmake with options: \"$2\""
		cmake "${REPODIR}/." "-DCMAKE_BUILD_TYPE=Release" "$2"
	fi

	# Compile the library
	echo "Building the GraphUtilities library"
	make -j

	# Run the tests
	echo "Running the validation tests"
	ctest

	# Check if the validation was successfull and exit if not
	if [[ $? != 0 ]]; then
		echo "Validation failed for build: $1"
		exit 1
	fi
}


# --------------------
# Sub-build n°1 : bare
# --------------------

runSubBuild "default" ""


# --------------------------
# Sub-build n°2 : STATISTICS
# --------------------------

runSubBuild "benchmarks" "-DENABLE_STATISTICS=ON"


# --------------------------
# Sub-build n°3 : BENCHMARKS
# --------------------------

runSubBuild "benchmarks" "-DENABLE_BENCHMARKS=ON"


# -------------------
# Sub-build n°4 : no-TLS
# -------------------

runSubBuild "tls" "-DENABLE_TLS=OFF"



# Exit with exit-code 0 if all sub-builds were successfull
exit 0
