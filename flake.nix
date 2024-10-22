# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
#
# SPDX-License-Identifier: AGPL-3.0-only
{
  outputs = {
    self,
    nixpkgs,
    flake-parts,
    pre-commit-hooks-nix,
    std,
    ...
  } @ inputs: let
    supportedSystems = [
      "x86_64-linux"
      "aarch64-linux"
      "x86_64-darwin"
      "aarch64-darwin"
    ];
    nixpkgsConfig = {
      allowUnsupportedSystem = true;
    };
    std_flake_out =
      std.growOn {
        inherit inputs nixpkgsConfig;
        cellsFrom = ./nix;

        cellBlocks = with std.blockTypes; [
          (runnables "apps")
          (installables "packages")
          (pkgs "pkgs")
          (devshells "devshells")
        ];
      }
      {
        packages = std.harvest inputs.self [["local" "packages"] ["local" "pkgs"]];
        devShells = std.harvest inputs.self ["local" "devshells"];
      };
    flake_parts_out = flake-parts.lib.mkFlake {inherit inputs;} {
      flake = {inherit std_flake_out;};
      imports = [
        inputs.pre-commit-hooks-nix.flakeModule
        inputs.std.flakeModule
        inputs.treefmt-nix.flakeModule
      ];
      systems = supportedSystems;
      perSystem = {
        config,
        inputs',
        lib,
        pkgs,
        system,
        ...
      }: rec {
        _module.args.pkgs = import nixpkgs {
          inherit system;
          config = nixpkgsConfig;
        };
        pre-commit = {
          check.enable = false; # Do not run in nix flake check because our dependencies are not available in this environment.
          settings = {
            # addGcRoot fixes issue:
            # =====> .pre-commit-config.yaml is not a file
            # introduced in
            # https://github.com/cachix/git-hooks.nix/commit/55b98216505b209b09499cfa39b34a3d63bc9ca3
            # to enforce line 367
            # ln -fs ${configFile} "''${GIT_WC}/.pre-commit-config.yaml"
            addGcRoot = false;
            hooks = {
              # C
              clang-format.enable = true;
              clang-tidy.enable = true;

              # Nix
              alejandra.enable = true;
              deadnix = {
                enable = true;
                settings = {
                  noUnderscore = true;
                  noLambdaPatternNames = true;
                };
              };
              statix.enable = true;

              # Python
              black.enable = true;
              isort.enable = true;
              mypy.enable = true;
              pylint.enable = true;
              pyupgrade.enable = true;

              # Misc
              check-added-large-files.enable = true;
              check-case-conflicts.enable = true;
              check-executables-have-shebangs.enable = true;
              check-shebang-scripts-are-executable.enable = true;
              check-symlinks.enable = true;
              check-vcs-permalinks.enable = true;
              detect-private-keys.enable = true;
              end-of-file-fixer.enable = true;
              mixed-line-endings.enable = true;
              reuse.enable = true;
              trim-trailing-whitespace.enable = true;
              typos = {
                enable = true;
                settings = {
                  configuration = ''
                    [default]
                    check-filename = false
                    extend-ignore-re = [
                      "PNGs",
                      "ba",
                      "ND",
                      'CompressedData\["(.|\n)*"\]',
                    ]
                  '';
                  locale = "en-us";
                };
              };

              # TOML
              check-toml.enable = true;
              taplo.enable = true;
            };
          };
        };
        devShells = let
          inherit (std_flake_out.packages.${system}) epsteinlib;
        in rec {
          epstein_devshell = pkgs.mkShell {
            inherit (config.pre-commit.devShell) shellHook nativeBuildInputs;

            inputsFrom = [
              std_flake_out.devShells.${system}.std
              epsteinlib
            ];
            packages = [
              epsteinlib
            ];
          };
          default = epstein_devshell;
        };
        packages = std_flake_out.packages.${system};

        treefmt = {
          projectRootFile = "flake.nix";
          programs = {
            # JSON
            biome.enable = true;
            # C
            clang-format.enable = true;
            # Nix
            alejandra.enable = true;
            # Python
            black.enable = true;
            # TOML
            taplo.enable = true;
            # Yaml
            yamlfmt.enable = true;
          };
        };
      };
    };
  in
    std_flake_out // flake_parts_out;

  inputs = {
    nixpkgs-stable.url = "github:nixos/nixpkgs/nixos-24.05";
    nixpkgs-unstable.url = "github:nixos/nixpkgs/nixos-unstable";
    nixpkgs.follows = "nixpkgs-stable";

    flake-parts.url = "github:hercules-ci/flake-parts";

    pre-commit-hooks-nix.url = "github:cachix/git-hooks.nix";

    treefmt-nix = {
      url = "github:numtide/treefmt-nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    devshell = {
      url = "github:numtide/devshell";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    std = {
      url = "github:divnix/std";
      inputs = {
        nixpkgs.follows = "nixpkgs";
        devshell.url = "github:numtide/devshell";
        devshell.follows = "devshell";
      };
    };
  };
}
