To run a test locally, first build the project

```bash
cmake -S . -B build
cmake --build build --parallel
```

Then enter the build folder
```bash
cd ./build
```

Then run `ctest`
```bash
ctest
```

And you're done!