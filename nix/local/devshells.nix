{
  inputs,
  cell,
}: let
  inherit (inputs) nixpkgs std;
  l = nixpkgs.lib // builtins;
in
  l.mapAttrs (_: std.lib.dev.mkShell) {
    std = {...}: {
      name = "epstein devshell";
      imports = [std.std.devshellProfiles.default];
      packages = with nixpkgs; [git meson ninja pkg-config lapack doxygen_gui];

      commands = [
        {
          name = "tests";
          command = "echo not implemented...";
          help = "run the unit tests";
          category = "Testing";
        }
        {
          name = "format";
          command = "nix fmt";
          help = "format all files";
          category = "Tooling";
        }
        {
          name = "docs";
          command = "doxygen Doxyfile && ${nixpkgs.xdg-utils}/bin/xdg-open html/index.html";
          help = "generate and show documentation";
          category = "Tooling";
        }
      ];
    };
  }
