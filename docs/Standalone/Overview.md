# Standalone Calculations

## Overview

In the workflow of designing a QM/MM umbrella sampling automation program, the tools required to run standard MD and QM/MM calculations were generated. 

We therefore added an interface that access these tools to allow the end user to quickly setup and run MD and QM/MM calculations. The tool uses the same file formats as the umbrella sampling code, and also reads in standard [MD parameters](../UserVars.md#md-global-variables) that can be edited in the [UserVars](../UserVars.md) directory.  

## Usage

To use the standalone code, you can run:

``` sh 
Standalone --help # Gets all cli variables
```

For a more in-depth usage guide: [Read this](./usage.md)

## Examples

There are also [example input files](./Examples.md) for a variety of different calculation types. 