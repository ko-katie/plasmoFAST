export function extractSequenceLines(chunks: string[]): string[] {
  const sequences: string[] = [];
  let leftover = '';
  let lineIndex = 0;

  for (const chunk of chunks) {
    const text = leftover + chunk;
    const lines = text.split('\n');
    leftover = lines.pop() ?? '';

    for (const line of lines) {
      if (line.trim() === '') continue;
      if (lineIndex % 4 === 1) sequences.push(line);
      lineIndex++;
    }
  }

  if (leftover.trim() !== '') {
    if (lineIndex % 4 === 1) sequences.push(leftover);
    lineIndex++;
  }

  return sequences;
}

export function isGzip(file: File): boolean {
  return file.name.endsWith('.gz');
}

export function getTextStream(file: File): ReadableStream<string> {
  let stream = file.stream() as unknown as ReadableStream<Uint8Array>;
  if (isGzip(file)) {
    stream = stream.pipeThrough(new DecompressionStream('gzip') as unknown as TransformStream<Uint8Array, Uint8Array>);
  }
  return stream.pipeThrough(new TextDecoderStream() as unknown as TransformStream<Uint8Array, string>);
}
