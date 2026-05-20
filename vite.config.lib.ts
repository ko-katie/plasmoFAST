import { defineConfig } from 'vite';
import { resolve } from 'path';
import { fileURLToPath } from 'url';
import dts from 'vite-plugin-dts';

const __dirname = fileURLToPath(new URL('.', import.meta.url));

export default defineConfig({
  plugins: [
    dts({ include: ['src'], exclude: ['src/**/*.test.ts'] }),
  ],
  build: {
    lib: {
      entry: resolve(__dirname, 'src/index.ts'),
      formats: ['es'],
      fileName: 'index',
    },
    rollupOptions: {
      input: {
        index: resolve(__dirname, 'src/index.ts'),
        'counter.worker': resolve(__dirname, 'src/counter.worker.ts'),
      },
      output: {
        entryFileNames: '[name].js',
        preserveModules: false,
      },
    },
  },
});
