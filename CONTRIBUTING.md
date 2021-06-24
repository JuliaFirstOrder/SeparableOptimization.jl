## Filing issues

When filing an issue, make sure to answer these five questions:

1. What version of Julia are you using (`julia --version`)?
2. What did you do?
3. What did you expect to see?
4. What did you see instead?

## Report a Bug

Open an issue. Please include descriptions of the following:
- The version of Julia that you are running
- The version of LCSO.jl that you are running
- Observations
- Expectations
- Steps to reproduce

## Contributing code

In general, this project follows Rust project conventions. Please make sure
you've linted, formatted, and run your tests before submitting a patch.

## Contribute a Bug Fix

- Report the bug first
- Create a pull request for the fix

## Suggest a New Feature

- Create a new issue to start a discussion around new topic. Label the issue as `new-feature`

# Developer guidelines

## Code formatting
All code should be autoformatted using [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl).  We use the options contained in `.JuliaFormatter.toml`. To format the code, `cd` to the project root, then:
1. install `] add JuliaFormatter`
2. format
  ```julia
  using JuliaFormatter
  format("src")
  format("test")
  ```

## Code coverage
Code coverage can be generated locally using the [LocalCoverage.jl](https://github.com/JuliaCI/LocalCoverage.jl) package.
1. install `] add LocalCoverage`
2. generate coverage
  ```julia
  using LocalCoverage
  generate_coverage("LCSO")
  ```
3. open coverage report html `open_coverage("LCSO")`
