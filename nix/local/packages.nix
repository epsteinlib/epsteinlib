# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
#
# SPDX-License-Identifier: AGPL-3.0-only
{
  inputs,
  cell,
}: let
  inherit (inputs) nixpkgs;
  pkgs = nixpkgs;
  pythonPackages = with pkgs.python3Packages; [cython numpy mpmath];
  build_epstein = buildtype: (pkgs.stdenv.mkDerivation {
    name = "libepstein";
    src = inputs.self;
    outputs = ["out" "dev"];

    nativeBuildInputs = with pkgs; [meson ninja pkg-config python3] ++ pythonPackages;
    buildInputs = [];
    mesonBuildType = buildtype;

    enableParallelBuilding = true;

    meta = {
      homepage = "https://feanor.num.uni-sb.de/gutendorf/epstein_bibliothek";
      # license = with licenses; [gpl3Only];
      mainProgram = "epsteinlib_c-lattice_sum";
    };
  });
in rec {
  epsteinlib = build_epstein "release";
  epsteinlib_dbg = build_epstein "debug";
  default = epsteinlib;
}
