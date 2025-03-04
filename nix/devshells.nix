# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
#
# SPDX-License-Identifier: AGPL-3.0-only
_: {
  perSystem = {pkgs, ...}: let
    enter_project_root_cmd = "pushd $(git rev-parse --show-toplevel)";
    docs_gen_cmd = "${enter_project_root_cmd} &&
                  mkdir -p html &&
                  PROJECT_NUMBER=$__VERSION doxygen Doxyfile &&
                  popd";
  in {
    devshells._epstein = {
      devshell.prj_root_fallback.eval = "$(git rev-parse --show-toplevel)";
      commands = [
        {
          name = "tests";
          command = "${enter_project_root_cmd} &&
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
          command = ''
            ${enter_project_root_cmd}/python &&
            stubgen epsteinlib.pyx -o .out &&
            mv .out/__main__.pyi epsteinlib.pyi &&
            rm -r .out &&
            format &&
            popd
          '';
          help = "regenerate the python stub files";
          category = "Tooling";
        }
        {
          name = "docs";
          command = ''
            ${enter_project_root_cmd} &&
            __VERSION=$(cat VERSION) &&
            ${docs_gen_cmd} &&
            ln -sf $__VERSION html/latest &&
            (${pkgs.xdg-utils}/bin/xdg-open html/latest/index.html || true)
            popd
          '';
          help = "generate and show documentation";
          category = "Tooling";
        }
        {
          name = "docs_all";
          command = ''
            docs
            for tag in $(git tag -l 'v*'); do
              git checkout $tag
              git checkout - Doxyfile # Use latest Doxyfile
              __VERSION=''${tag:1}
              ${docs_gen_cmd}
              git checkout - # Back to original checkout
            done
          '';
          help = "generate latest and documentation for previous releases";
          category = "Tooling";
        }
        {
          name = "website";
          command = ''
            ${enter_project_root_cmd}/website &&
            ln -sf ../html .
            hugo server --buildDrafts --disableFastRender --renderToMemory
            popd
          '';
          help = "Show website with live editing";
          category = "Tooling";
        }
        {
          name = "website_deploy";
          command = ''
            ${enter_project_root_cmd}/website &&
            ln -sf ../html .
            HUGO_ENVIRONMENT=production hugo --gc --minify &&
            find public -type f -name "*.license" -delete
            GH_PAGES_BRANCH=gh-pages
            git worktree add $GH_PAGES_BRANCH || true
            rm -rf $GH_PAGES_BRANCH/* && mv public/* $GH_PAGES_BRANCH/ &&
            pushd $GH_PAGES_BRANCH
            git add . && git commit -m 'Deploy website' &&
            git push origin $GH_PAGES_BRANCH &&
            popd
            git worktree remove $GH_PAGES_BRANCH
            popd
          '';
          help = "Builds static website in website/public folder and pushes it to gh-pages branch for deployment";
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
