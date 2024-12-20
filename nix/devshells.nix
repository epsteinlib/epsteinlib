# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
#
# SPDX-License-Identifier: AGPL-3.0-only
_: {
  perSystem = {
    pkgs,
    self',
    ...
  }: {
    devshells._epstein = {
      devshell.prj_root_fallback.eval = "$(git rev-parse --show-toplevel)";
      commands = [
        {
          name = "tests";
          command = "pushd $(git rev-parse --show-toplevel) &&
                     meson setup --reconfigure build -Db_coverage=true &&
                     meson compile -C build &&
                     meson test -v -C build $@
                     mkdir -p html &&
                     gcovr --html-details html/coverage.html --txt --txt-metric branch --print-summary --exclude 'build/python/.*pyx.c' --exclude 'src/tests/.*.c' &&
                     popd
                    ";
          help = "run the unit tests (arguments are passed to meson test, e.g. for only running specific tests)";
          category = "Testing";
        }
        {
          name = "format";
          command = "nix fmt";
          help = "format all files";
          category = "Tooling";
        }
        {
          name = "generate_python_stubs";
          command = "pushd $(git rev-parse --show-toplevel)/python &&
                     stubgen epsteinlib.pyx -o .out &&
                     mv .out/__main__.pyi epsteinlib.pyi &&
                     rm -r .out &&
                     nix fmt &&
                     popd";
          help = "regenerate the python stub files";
          category = "Tooling";
        }
        {
          name = "docs";
          command = "PROJECT_NUMBER=$(cat VERSION) doxygen Doxyfile && (${pkgs.xdg-utils}/bin/xdg-open html/index.html || true)";
          help = "generate and show documentation";
          category = "Tooling";
        }
        {
          name = "build_release";
          command = "rm -rf dist && pyproject-build --sdist";
          help = "Build release which creates dists folder";
          category = "Releasing";
        }
        {
          name = "testpypi_upload";
          command = "twine upload --repository testpypi dist/*";
          help = "Upload dist/* folder to https://tests.pypi.org";
          category = "Releasing";
        }
        {
          name = "testpypi_install";
          command = "python -m pip install --force -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ epsteinlib";
          help = "Install latest uploaded version from https://tests.pypi.org";
          category = "Releasing";
        }
      ];
    };
  };
}
