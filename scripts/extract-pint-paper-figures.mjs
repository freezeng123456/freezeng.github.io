import { execFileSync } from "node:child_process"
import { mkdirSync, readFileSync, unlinkSync, writeFileSync } from "node:fs"
import { basename, dirname, join, resolve } from "node:path"

const paperPath = process.argv[2]
const pdftocairo = process.argv[3] ?? process.env.PDFTOCAIRO ?? "pdftocairo"

if (!paperPath) {
  console.error("Usage: node scripts/extract-pint-paper-figures.mjs <paper.pdf> [pdftocairo]")
  process.exit(1)
}

const outputDirectory = resolve("content/assets/papers/time-parallelization/source-figures")
mkdirSync(outputDirectory, { recursive: true })

// Coordinates are in PDF points. Most pages are 493.228 x 700.157 points;
// physical PDF page 95 (Table 4.1) is a 700.157 x 493.228 landscape page. Each
// crop keeps the complete graphical or tabular object and omits the running
// header and body prose. Captions are translated next to the asset in the
// bilingual notes.
const assets = [
  ["figure-1-1", 3, [44, 54, 405, 117]],
  ["figure-2-1", 6, [44, 43, 405, 228]],
  ["figure-2-2", 8, [44, 43, 405, 440]],
  ["figure-2-3", 10, [44, 43, 405, 440]],
  ["figure-2-4", 11, [44, 43, 405, 230]],
  ["figure-3-1", 17, [44, 43, 405, 175]],
  ["figure-3-2", 19, [44, 43, 405, 308]],
  ["figure-3-3", 21, [44, 43, 405, 220]],
  ["figure-3-4", 24, [44, 43, 405, 362]],
  ["figure-3-5", 26, [44, 43, 405, 346]],
  ["figure-3-6", 27, [44, 43, 405, 302]],
  ["figure-3-7", 28, [44, 43, 405, 137]],
  ["figure-3-8", 31, [44, 43, 405, 121]],
  ["figure-3-9", 36, [44, 43, 405, 160]],
  ["figure-3-10", 37, [44, 43, 405, 160]],
  ["figure-3-11", 39, [44, 43, 405, 178]],
  ["table-3-1", 40, [44, 50, 405, 107]],
  ["figure-3-12", 43, [44, 43, 405, 175]],
  ["figure-3-13", 44, [44, 43, 405, 170]],
  ["table-3-2", 45, [44, 50, 405, 178]],
  ["figure-3-14", 46, [44, 43, 405, 170]],
  ["figure-3-15", 50, [44, 43, 405, 486]],
  ["figure-3-16", 53, [44, 43, 405, 160]],
  ["figure-3-17", 57, [44, 43, 405, 249]],
  ["figure-3-18", 58, [44, 43, 405, 160]],
  ["figure-4-1", 61, [44, 43, 405, 90]],
  ["figure-4-2", 64, [44, 43, 405, 160]],
  ["figure-4-3", 66, [44, 43, 405, 120]],
  ["figure-4-4", 67, [44, 43, 405, 113]],
  ["figure-4-5", 67, [44, 205, 405, 172]],
  ["figure-4-6", 71, [44, 43, 405, 233]],
  ["figure-4-7", 72, [44, 43, 405, 82]],
  ["figure-4-8", 73, [44, 43, 405, 259]],
  ["figure-4-9", 75, [44, 43, 405, 245]],
  ["figure-4-10", 76, [44, 43, 405, 307]],
  ["figure-4-11", 77, [44, 43, 405, 118]],
  ["figure-4-12", 80, [44, 43, 405, 178]],
  ["figure-4-13", 81, [44, 43, 405, 162]],
  ["figure-4-14", 85, [44, 43, 405, 162]],
  ["figure-4-15", 86, [44, 43, 405, 178]],
  ["figure-4-16", 86, [44, 249, 405, 163]],
  ["figure-4-17", 87, [44, 43, 405, 178]],
  ["figure-4-18", 92, [44, 43, 405, 126]],
  ["figure-4-19", 93, [44, 43, 405, 178]],
  ["figure-4-20", 93, [44, 274, 405, 179]],
  ["figure-4-21", 94, [44, 43, 405, 251]],
  ["table-4-1", 95, [84, 54, 562, 334]],
  ["figure-4-22", 97, [44, 43, 405, 173]],
]

function escapeXml(value) {
  return value
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&apos;")
}

for (const [name, pdfPage, [x, y, width, height]] of assets) {
  const rawPath = join(outputDirectory, `.${name}-page.svg`)
  const outputPath = join(outputDirectory, `${name}.svg`)

  execFileSync(
    pdftocairo,
    ["-svg", "-f", String(pdfPage), "-l", String(pdfPage), paperPath, rawPath],
    { stdio: "inherit" },
  )

  const articlePage = pdfPage + 384
  const printedLabel = name.startsWith("table")
    ? `Table ${name.slice("table-".length).replace("-", ".")}`
    : `Figure ${name.slice("figure-".length).replace("-", ".")}`
  const titleId = `${name}-title`
  const descId = `${name}-desc`
  const title = `${printedLabel} from the source paper`
  const description = `${printedLabel} as printed on article page ${articlePage} of Time parallelization for hyperbolic and parabolic problems.`

  let svg = readFileSync(rawPath, "utf8")
  svg = svg.replace(/<svg\b([^>]*)>/, (_, attributes) => {
    const preserved = attributes
      .replace(/\s(?:width|height|viewBox|role|aria-labelledby|overflow)="[^"]*"/g, "")
      .trim()
    return `<svg ${preserved} width="${width}pt" height="${height}pt" viewBox="${x} ${y} ${width} ${height}" role="img" aria-labelledby="${titleId} ${descId}" overflow="hidden">`
  })
  svg = svg.replace(
    /(<svg\b[^>]*>)/,
    `$1\n<title id="${titleId}">${escapeXml(title)}</title>\n<desc id="${descId}">${escapeXml(description)}</desc>\n<rect x="${x}" y="${y}" width="${width}" height="${height}" fill="#ffffff"/>`,
  )
  svg = svg.replace(
    /<\?xml version="1\.0" encoding="UTF-8"\?>/,
    '<?xml version="1.0" encoding="UTF-8"?>\n<!-- Extracted from the CC BY 4.0 source paper; see the adjacent note attribution. -->',
  )

  writeFileSync(outputPath, svg)
  unlinkSync(rawPath)
  console.log(`${basename(outputPath)} <- PDF page ${pdfPage}`)
}

console.log(`Wrote ${assets.length} assets to ${dirname(join(outputDirectory, "x"))}`)
