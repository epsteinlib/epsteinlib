# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
#
# SPDX-License-Identifier: AGPL-3.0-only
{self, ...}: {
  perSystem = {
    pkgs,
    lib,
    self',
    ...
  }: let
    pythonDevEnv = pkgs.python3.withPackages (
      ps:
        with ps; [
          twine
          build
          mypy
          mpmath
          numpy
          matplotlib
          self'.packages.epsteinlib_python
          pylint
        ]
    );

    build_epstein = buildtype: (pkgs.stdenv.mkDerivation {
      name = "epsteinlib";
      src = self;
      outputs = [
        "out"
        "dev"
      ];

      nativeBuildInputs = with pkgs; [
        meson
        ninja
        pkg-config
        (pkgs.python3.withPackages (
          ps:
            with ps; [
              cython
            ]
        ))
      ];
      buildInputs = [];
      mesonBuildType = buildtype;
      mesonFlags = ["-Dbuild_python=false"];

      enableParallelBuilding = true;

      passthru.optional-dependencies.dev = with pkgs; [
        # Linters/formatters
        git
        doxygen_gui
        graphviz
        neovim
        gcovr
        pythonDevEnv
        clang-tools
      ];
      meta = {
        homepage = "https://github.com/epsteinlib/epsteinlib";
        license = with lib.licenses; [agpl3Only];
        mainProgram = "epsteinlib_c-lattice_sum";
      };
    });
    build_epstein_python = pkgs.python3Packages.buildPythonPackage rec {
      name = "epsteinlib";
      src = self;
      outputs = [
        "out"
        "dev"
      ];

      build-system = with pkgs; [
        python3Packages.meson-python
        python3Packages.cython
        pkg-config
      ];
      dependencies = with pkgs; [
        python3Packages.numpy
      ];
      pyproject = true;
      buildInputs = [];

      nativeCheckInputs = with pkgs; [
        python3Packages.unittestCheckHook
        python3Packages.mpmath
      ];

      unittestFlagsArray = [
        "-s"
        "python/tests"
        "-v"
      ];

      enableParallelBuilding = true;
      pythonImportsCheck = [name];

      meta = {
        homepage = "https://github.com/epsteinlib/epsteinlib";
        license = with lib.licenses; [agpl3Only];
      };
    };
  in {
    packages = rec {
      epsteinlib = build_epstein "release";
      epsteinlib_dbg = build_epstein "debug";
      epsteinlib_python = build_epstein_python;
      default = epsteinlib;
      inherit pythonDevEnv;
    };
  };
}
