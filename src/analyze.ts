import type { AnalysisResult, AnalyzeOptions } from './types';

const VALID_EXTENSIONS = ['.fastq', '.fastq.gz', '.fq', '.fq.gz'];

export function analyze(file: File, options: AnalyzeOptions = {}): Promise<AnalysisResult> {
  const { referenceUrl, onProgress } = options;

  if (!VALID_EXTENSIONS.some(ext => file.name.endsWith(ext))) {
    return Promise.reject(
      new Error(`Unsupported file type: "${file.name}". Expected .fastq or .fastq.gz`)
    );
  }

  return new Promise((resolve, reject) => {
    const worker = new Worker(
      /* @vite-ignore */
      new URL('./counter.worker.js', import.meta.url),
      { type: 'module' }
    );

    worker.onmessage = (event: MessageEvent) => {
      const { type, data, message } = event.data;
      if (type === 'progress') {
        onProgress?.({ bytesRead: event.data.bytesRead, totalBytes: event.data.totalBytes });
      } else if (type === 'result') {
        worker.terminate();
        resolve(data as AnalysisResult);
      } else if (type === 'error') {
        worker.terminate();
        reject(new Error(message));
      }
    };

    worker.onerror = (event: ErrorEvent) => {
      worker.terminate();
      reject(event);
    };

    worker.postMessage({ file, referenceUrl });
  });
}
