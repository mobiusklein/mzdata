import { test_reader } from "mzdata-wasm";

async function readFile(handle: File) {
  const stream = await handle.stream();
  const reader = stream.getReader();
  return await test_reader(reader);
}

export { readFile };
