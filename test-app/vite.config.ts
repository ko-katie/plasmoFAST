import { defineConfig } from 'vite';
import { resolve } from 'path';
import { fileURLToPath } from 'url';

const __dirname = fileURLToPath(new URL('.', import.meta.url));

export default defineConfig({
  root: __dirname,
  resolve: {
    alias: {
      'plasmofast': resolve(__dirname, '../src/index.ts'),
    },
  },
});
