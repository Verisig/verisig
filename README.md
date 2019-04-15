# Verisig: a tool for verifying safety properties of closed-loop systems with neural network components
## Verisig Source Code: https://github.com/verisig/verisig

***********************************
### Contributors
***********************************
* Radoslav Ivanov
* Taylor J. Carpenter
* James Weimer
* Rajeev Alur
* George J. Pappas
* Insup Lee

***********************************
### References
***********************************
The primary Verisig paper is this one:
* R. Ivanov, J. Weimer, R. Alur, G. J. Pappas, and I. Lee. Verisig: verifying safety properties of hybrid systems with neural network controllers. _Proceedings of the 22nd International Converence on Hybrid Systems: Computation and Control_, 2019 [https://arxiv.org/abs/1811.01828]
***********************************
### VERISIG USAGE
***********************************
Verisig has been tested on Ubuntu (16.04 and 18.04) using Java 1.8. A complete explanation of the Verisig workflow can be found in the user manual.

After installing the distribution using `./gradlew installDist`, you can use the executor script with the -h flag to see available command-line options:
```
$ ./verisig -h
Verisig v0.9 Usage:
verisig [options...] spaceex_xml dnn_yml
 --flowstar-cmd CMD             : Custom flow* command (default: flowstar)
 --help (-h)                    : print help message (default: true)
 --no-flowstar (-nf)            : System composition only (default: false)
 --output-model (-o)            : Enable flow* model output (default: false)
 --output-model-name (-of) PATH : Flow* model output location (overrides -o)
 --spaceex-config (-sc) PATH    : SpaceEx 'cfg' file
 --verbose (-v)                 : enable verbose printing (default: false)
 --verisig-config (-vc) PATH    : Verisig yaml config
 --version                      : print version (default: false)
```

`--flowstar-cmd CMD` is used to specify the path to the Flow* executable. The default value is `flowstar`, meaning the executable is on the PATH. This value is not used if `--no-flowstar` is present.

`--no-flowstar` is used to indicate that Flow* should not be run. When this flag is present, Verisig composes the hybrid system model but does not verify the system.

`--output-model` is used to indicate that the composed model file should be written to disk. Without this flag, the model is only kept in memory. The model name matches that of the SpaceEx file name but with a _.model_ extension. When combined with the `--no-flowstar` flag, Verisig writes out the model to be inspected or run later.

`--output-model-name PATH` is used to specify the path to use for output the composed Flow* model. This option forces the `--output-model` flag.

`--spaceex-config PATH` is used to specify the path to the SpaceEx _.cfg_ file. The default value is that of the SpaceEx XML file but with a _.cfg_ extension.

`--verisig-config PATH` is used to specify the path to the Verisig _.yml_ configuration file. The default value is that of the SpaceEx XML file but with a _.yml_ extension.


#### Example
The following command can be run from the local directory to verify a single interval of the cartpole example:

```
./verisig --flowstar-command ./flowstar/flowstar examples/mountain_car/MC.xml examples/mountain_car/sig16x16.yml
```

Look at the file __examples/mountain_car/multi_runner.py__ to see how one can verify the entire range of the unsafe set that was used in the case-study.