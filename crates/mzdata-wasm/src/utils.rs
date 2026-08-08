

// const subsampleResolutionSpacing = (x: number[], desiredResolution: number) => {
//     const keptIndices = [0];
//     if (x.length == 0) return keptIndices

//     let last = x[0]
//     for (let i = 1; i < x.length; i++) {
//         if (x[i] - last > desiredResolution) {
//             keptIndices.push(i);
//             last = x[i]
//         }
//     }
//     if (keptIndices[keptIndices.length - 1] != x.length - 1) {
//         keptIndices.push(x.length - 1);
//     }
//     return keptIndices;
// };
