# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
#
# SPDX-License-Identifier: AGPL-3.0-only
{
  outputs = {
    self,
    nixpkgs,
    flake-parts,
    pre-commit-hooks-nix,
    devshell,
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
    flake_parts_out = flake-parts.lib.mkFlake {inherit inputs;} {
      imports = [
        inputs.pre-commit-hooks-nix.flakeModule
        inputs.treefmt-nix.flakeModule
        inputs.devshell.flakeModule
        ./nix/packages.nix
        ./nix/devshells.nix
      ];
      systems = supportedSystems;
      perSystem = {
        self',
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

            hooks = let
              mypy_wrapper = pkgs.writeShellScript "mypy" ''
                MYPYPATH=${self'.packages.epsteinlib_python}/lib/python3.11/site-packages/ ${self'.packages.pythonDevEnv}/bin/mypy "$@"
              '';
              pylint_wrapper = pkgs.writeShellScript "pylint" ''
                ${self'.packages.pythonDevEnv}/bin/pylint "$@"
              '';
            in {
              # C
              clang-format.enable = true;
              clang-tidy = {
                enable = true;
              };

              # Nix
              alejandra.enable = false;
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
              mypy = {
                enable = true;
                settings.binPath = "${mypy_wrapper}";
              };
              pylint = {
                enable = true;
                settings.binPath = "${pylint_wrapper}";
              };
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
          inherit (self'.devShells) _epstein;
        in rec {
          epstein = pkgs.mkShell {
            inherit (_epstein) shellHook;

            inputsFrom = [
              self'.packages.default
            ];
            packages = [
              self'.packages.default
              self'.packages.default.optional-dependencies.dev
            ];
          };
          epstein_devshell = pkgs.mkShell {
            shellHook = config.pre-commit.installationScript + epstein.shellHook;

            inputsFrom = [
              epstein
            ];
            packages = [
              config.pre-commit.settings.enabledPackages
            ];
          };
          default = epstein_devshell;
        };

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
    flake_parts_out;

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
  };
}
