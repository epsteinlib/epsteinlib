{
  inputs,
  cell,
}: let
  inherit (inputs) nixpkgs;
  pkgs = nixpkgs;
  build_epstein = buildtype: (pkgs.stdenv.mkDerivation {
    name = "libepstein";
    src = inputs.self;
    outputs = ["out" "dev"];

    nativeBuildInputs = with pkgs; [meson ninja pkg-config];
    buildInputs = [];
    mesonBuildType = buildtype;

    enableParallelBuilding = true;

    meta = {
      homepage = "https://feanor.num.uni-sb.de/gutendorf/epstein_bibliothek";
      # license = with licenses; [gpl3Only];
      mainProgram = "test_exec";
    };
  });
in rec {
  epsteinlib = build_epstein "release";
  epsteinlib_dbg = build_epstein "debug";
  default = epsteinlib;
}
