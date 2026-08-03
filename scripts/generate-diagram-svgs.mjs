import fs from "node:fs"
import path from "node:path"
import { fileURLToPath } from "node:url"

const scriptDir = path.dirname(fileURLToPath(import.meta.url))
const repoRoot = path.resolve(scriptDir, "..")
const assetRoot = path.join(repoRoot, "content", "assets", "diagrams")

const C = {
  ink: "#14213d",
  muted: "#5e6b7f",
  paper: "#fffdf8",
  line: "#b8c4d2",
  teal: "#0f766e",
  tealSoft: "#dff6f1",
  blue: "#2563eb",
  blueSoft: "#e9f0ff",
  indigo: "#4f46e5",
  indigoSoft: "#eeecff",
  amber: "#d97706",
  amberSoft: "#fff2d9",
  rose: "#dc4c64",
  roseSoft: "#ffe9ed",
  green: "#16865b",
  greenSoft: "#e2f7ec",
  slateSoft: "#edf2f7",
  white: "#ffffff",
}

const esc = (value) =>
  String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")

function textBlock(x, y, lines, options = {}) {
  const {
    size = 22,
    weight = 600,
    fill = C.ink,
    anchor = "middle",
    lineHeight = Math.round(size * 1.35),
    family = "Inter, Noto Sans SC, PingFang SC, Microsoft YaHei, system-ui, sans-serif",
    opacity = 1,
  } = options
  const start = y - ((lines.length - 1) * lineHeight) / 2
  return `<text x="${x}" y="${start}" text-anchor="${anchor}" font-family="${family}" font-size="${size}" font-weight="${weight}" fill="${fill}" opacity="${opacity}">${lines
    .map(
      (line, index) => `<tspan x="${x}" dy="${index === 0 ? 0 : lineHeight}">${esc(line)}</tspan>`,
    )
    .join("")}</text>`
}

function pill(x, y, w, label, color = C.teal, fill = C.tealSoft, options = {}) {
  const { h = 34, size = 15 } = options
  return `<g>
    <rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${h / 2}" fill="${fill}" stroke="${color}" stroke-opacity=".28"/>
    ${textBlock(x + w / 2, y + h / 2 + size * 0.34, [label], { size, weight: 700, fill: color })}
  </g>`
}

function card(x, y, w, h, options = {}) {
  const {
    title,
    body = [],
    accent = C.teal,
    fill = C.white,
    step,
    titleSize = 21,
    bodySize = 16,
    align = "left",
    dashed = false,
  } = options
  const tx = align === "left" ? x + 28 : x + w / 2
  const anchor = align === "left" ? "start" : "middle"
  const titleY = body.length ? y + 43 : y + h / 2 + titleSize * 0.34
  const bodyY = y + 75
  return `<g filter="url(#shadow)">
    <rect x="${x}" y="${y}" width="${w}" height="${h}" rx="22" fill="${fill}" stroke="${accent}" stroke-opacity=".28" stroke-width="1.5" ${dashed ? 'stroke-dasharray="7 7"' : ""}/>
    <rect x="${x}" y="${y}" width="8" height="${h}" rx="4" fill="${accent}"/>
    ${
      step
        ? `<circle cx="${x + w - 28}" cy="${y + 28}" r="17" fill="${accent}"/>
           ${textBlock(x + w - 28, y + 33, [step], { size: 14, weight: 800, fill: C.white })}`
        : ""
    }
    ${textBlock(tx, titleY, Array.isArray(title) ? title : [title], { size: titleSize, weight: 750, fill: C.ink, anchor })}
    ${body.length ? textBlock(tx, bodyY, body, { size: bodySize, weight: 500, fill: C.muted, anchor, lineHeight: Math.round(bodySize * 1.45) }) : ""}
  </g>`
}

function arrowMarker(color, fallback = "arrow") {
  return (
    new Map([
      [C.teal, "arrowTeal"],
      [C.blue, "arrowBlue"],
      [C.indigo, "arrowIndigo"],
      [C.amber, "arrowAmber"],
      [C.rose, "arrowRose"],
      [C.green, "arrowGreen"],
      [C.line, "arrow"],
    ]).get(color) ?? fallback
  )
}

function pathArrow(d, options = {}) {
  const { color = C.line, width = 2.6, marker = "arrow", dashed = false, opacity = 1 } = options
  const lineWidth = Math.min(width, 3.2)
  const dash = dashed ? 'stroke-dasharray="6 9"' : ""
  const resolvedMarker = arrowMarker(color, marker)
  return `<g class="flow-connector" opacity="${opacity}">
    <path d="${d}" fill="none" stroke="${C.white}" stroke-opacity=".88" stroke-width="${lineWidth + 4.8}" stroke-linecap="round" stroke-linejoin="round" ${dash}/>
    <path d="${d}" fill="none" stroke="${color}" stroke-width="${lineWidth}" stroke-linecap="round" stroke-linejoin="round" marker-end="url(#${resolvedMarker})" ${dash}/>
  </g>`
}

function lineArrow(x1, y1, x2, y2, options = {}) {
  return pathArrow(`M ${x1} ${y1} L ${x2} ${y2}`, options)
}

function section(x, y, w, h, title, color, options = {}) {
  const { subtitle = "", fill = C.white } = options
  return `<g>
    <rect x="${x}" y="${y}" width="${w}" height="${h}" rx="26" fill="${fill}" stroke="${color}" stroke-opacity=".22" stroke-width="1.5"/>
    <rect x="${x}" y="${y}" width="${w}" height="52" rx="26" fill="${color}" fill-opacity=".1"/>
    <rect x="${x}" y="${y + 26}" width="${w}" height="26" fill="${color}" fill-opacity=".1"/>
    ${textBlock(x + 24, y + 33, [title], { size: 17, weight: 800, fill: color, anchor: "start" })}
    ${subtitle ? textBlock(x + w - 24, y + 33, [subtitle], { size: 13, weight: 600, fill: C.muted, anchor: "end" }) : ""}
  </g>`
}

function frame({ width, height, kicker, title, subtitle, body }) {
  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="img" aria-labelledby="diagram-title diagram-desc">
  <title id="diagram-title">${esc(title)}</title>
  <desc id="diagram-desc">${esc(subtitle)}</desc>
  <defs>
    <linearGradient id="paper" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0" stop-color="#fffdf8"/>
      <stop offset=".55" stop-color="#f9fbff"/>
      <stop offset="1" stop-color="#f3faf8"/>
    </linearGradient>
    <pattern id="dots" width="28" height="28" patternUnits="userSpaceOnUse">
      <circle cx="2" cy="2" r="1.2" fill="#9fb0c2" opacity=".18"/>
    </pattern>
    <filter id="shadow" x="-20%" y="-20%" width="140%" height="150%">
      <feDropShadow dx="0" dy="8" stdDeviation="9" flood-color="#183153" flood-opacity=".10"/>
    </filter>
    <marker id="arrow" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="userSpaceOnUse">
      <path d="M 2 1.5 L 10 6 L 2 10.5" fill="none" stroke="${C.line}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
    <marker id="arrowTeal" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="userSpaceOnUse">
      <path d="M 2 1.5 L 10 6 L 2 10.5" fill="none" stroke="${C.teal}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
    <marker id="arrowBlue" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="userSpaceOnUse">
      <path d="M 2 1.5 L 10 6 L 2 10.5" fill="none" stroke="${C.blue}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
    <marker id="arrowIndigo" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="userSpaceOnUse">
      <path d="M 2 1.5 L 10 6 L 2 10.5" fill="none" stroke="${C.indigo}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
    <marker id="arrowAmber" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="userSpaceOnUse">
      <path d="M 2 1.5 L 10 6 L 2 10.5" fill="none" stroke="${C.amber}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
    <marker id="arrowRose" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="userSpaceOnUse">
      <path d="M 2 1.5 L 10 6 L 2 10.5" fill="none" stroke="${C.rose}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
    <marker id="arrowGreen" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="userSpaceOnUse">
      <path d="M 2 1.5 L 10 6 L 2 10.5" fill="none" stroke="${C.green}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
  </defs>
  <rect x="8" y="8" width="${width - 16}" height="${height - 16}" rx="34" fill="url(#paper)" stroke="#d9e2ec" stroke-width="1.5"/>
  <rect x="8" y="8" width="${width - 16}" height="${height - 16}" rx="34" fill="url(#dots)"/>
  ${pill(54, 42, Math.max(118, kicker.length * 11 + 36), kicker, C.teal, C.tealSoft, { h: 32, size: 14 })}
  ${textBlock(54, 112, [title], { size: 32, weight: 800, fill: C.ink, anchor: "start" })}
  ${textBlock(54, 146, [subtitle], { size: 16, weight: 500, fill: C.muted, anchor: "start" })}
  ${body}
</svg>
`
}

function servingLoop(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "LLM 推理",
          title: "Token 级服务循环",
          subtitle: "调度、执行、采样与 KV 内存管理在每个 step 中闭环",
          scheduler: "Token 调度器",
          schedulerBody: ["请求状态", "Token 预算"],
          model: "展平模型执行",
          modelBody: ["Embedding", "N × Transformer"],
          sampler: "LM Head + 采样器",
          samplerBody: ["Gather hidden states", "生成 next-token IDs"],
          output: "完成输出",
          kv: "KV Block Manager",
          kvBody: ["分配 · 映射 · 回收"],
          no: "未完成：进入下一 step",
          yes: "完成",
        }
      : {
          kicker: "LLM INFERENCE",
          title: "The Token-Level Serving Loop",
          subtitle:
            "Scheduling, execution, sampling, and KV memory form one step-wise control loop",
          scheduler: "Token Scheduler",
          schedulerBody: ["Request state", "Token budget"],
          model: "Flattened Execution",
          modelBody: ["Embedding", "N × Transformer"],
          sampler: "LM Head + Sampler",
          samplerBody: ["Gather hidden states", "Produce next-token IDs"],
          output: "Completed Output",
          kv: "KV Block Manager",
          kvBody: ["Allocate · map · reclaim"],
          no: "unfinished: next step",
          yes: "finished",
        }
  const body = `
    ${card(70, 235, 270, 145, { title: t.scheduler, body: t.schedulerBody, accent: C.blue, fill: C.blueSoft, step: "1" })}
    ${card(465, 205, 280, 145, { title: t.model, body: t.modelBody, accent: C.indigo, fill: C.indigoSoft, step: "2" })}
    ${card(860, 235, 270, 145, { title: t.sampler, body: t.samplerBody, accent: C.teal, fill: C.tealSoft, step: "3" })}
    ${card(860, 470, 270, 90, { title: t.output, accent: C.green, fill: C.greenSoft, align: "center" })}
    ${card(465, 430, 280, 115, { title: t.kv, body: t.kvBody, accent: C.amber, fill: C.amberSoft, align: "center" })}
    ${lineArrow(340, 307, 455, 277, { color: C.blue, marker: "arrowBlue" })}
    ${lineArrow(745, 277, 850, 307, { color: C.indigo, marker: "arrowBlue" })}
    ${lineArrow(995, 380, 995, 458, { color: C.green, marker: "arrowTeal" })}
    ${pill(1012, 403, 84, t.yes, C.green, C.greenSoft, { h: 28, size: 12 })}
    ${pathArrow("M 885 258 C 770 120, 300 120, 205 225", { color: C.rose, marker: "arrowRose", width: 3 })}
    ${pill(470, 158, lang === "zh" ? 210 : 178, t.no, C.rose, C.roseSoft, { h: 30, size: 12 })}
    ${pathArrow("M 530 430 C 460 405, 390 365, 335 345", { color: C.amber, marker: "arrowAmber", dashed: true })}
    ${pathArrow("M 680 430 C 745 395, 785 335, 855 320", { color: C.amber, marker: "arrowAmber", dashed: true })}
    ${pill(523, 568, 164, "request state + KV", C.amber, C.amberSoft, { h: 30, size: 12 })}
  `
  return frame({
    width: 1200,
    height: 630,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function pagedKv(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "PAGEDATTENTION",
          title: "从逻辑 Token 到物理 KV Slot",
          subtitle: "Block table 把动态序列映射到可复用的物理缓存块",
          logical: "逻辑位置",
          virtual: "虚拟块 1",
          offset: "块内偏移 9",
          table: "Block Table",
          request: "请求",
          logicalBlock: "逻辑块",
          physicalBlock: "物理块",
          cache: "物理 KV Cache",
          block: "Block 8",
          formula: "物理 slot = 8 × 16 + 9 = 137",
          use: "在每一层模型的 K/V cache 中使用 slot 137",
        }
      : {
          kicker: "PAGEDATTENTION",
          title: "From a Logical Token to a Physical KV Slot",
          subtitle: "A block table maps a dynamic sequence onto reusable physical cache blocks",
          logical: "Logical position",
          virtual: "Virtual block 1",
          offset: "Offset 9",
          table: "Block Table",
          request: "Request",
          logicalBlock: "Logical block",
          physicalBlock: "Physical block",
          cache: "Physical KV Cache",
          block: "Block 8",
          formula: "physical slot = 8 × 16 + 9 = 137",
          use: "Use slot 137 in every model layer's K/V cache",
        }

  let virtualCells = ""
  for (let i = 0; i < 16; i++) {
    const col = i % 8
    const row = Math.floor(i / 8)
    const x = 295 + col * 27
    const y = 285 + row * 34
    const active = i === 9
    virtualCells += `<rect x="${x}" y="${y}" width="22" height="26" rx="6" fill="${active ? C.amber : C.white}" stroke="${active ? C.amber : C.line}" stroke-width="1.4"/>
      ${textBlock(x + 11, y + 18, [i], { size: 10, weight: 700, fill: active ? C.white : C.muted })}`
  }

  let cacheBlocks = ""
  for (let i = 0; i < 12; i++) {
    const col = i % 6
    const row = Math.floor(i / 6)
    const x = 870 + col * 38
    const y = 274 + row * 58
    const active = i === 8
    cacheBlocks += `<g>
      <rect x="${x}" y="${y}" width="31" height="42" rx="7" fill="${active ? C.teal : C.white}" stroke="${active ? C.teal : C.line}" stroke-width="1.4"/>
      ${textBlock(x + 15.5, y + 26, [i], { size: 11, weight: 750, fill: active ? C.white : C.muted })}
    </g>`
  }

  const body = `
    <g filter="url(#shadow)">
      <circle cx="120" cy="310" r="72" fill="${C.blueSoft}" stroke="${C.blue}" stroke-width="2"/>
      ${textBlock(120, 284, [t.logical], { size: 16, weight: 700, fill: C.blue })}
      ${textBlock(120, 336, ["25"], { size: 42, weight: 850, fill: C.ink })}
    </g>
    ${section(255, 200, 275, 210, t.virtual, C.amber, { subtitle: t.offset, fill: C.amberSoft })}
    ${virtualCells}
    ${section(570, 200, 235, 210, t.table, C.indigo, { fill: C.indigoSoft })}
    ${textBlock(594, 281, [t.request, t.logicalBlock, t.physicalBlock], { size: 14, weight: 650, fill: C.muted, anchor: "start", lineHeight: 39 })}
    ${textBlock(775, 281, ["r", "1", "8"], { size: 15, weight: 800, fill: C.indigo, anchor: "end", lineHeight: 39, family: "IBM Plex Mono, ui-monospace, monospace" })}
    <line x1="590" y1="300" x2="785" y2="300" stroke="${C.indigo}" stroke-opacity=".18"/>
    <line x1="590" y1="339" x2="785" y2="339" stroke="${C.indigo}" stroke-opacity=".18"/>
    ${section(845, 200, 300, 210, t.cache, C.teal, { subtitle: t.block, fill: C.tealSoft })}
    ${cacheBlocks}
    ${lineArrow(194, 310, 244, 310, { color: C.blue, marker: "arrowBlue" })}
    ${lineArrow(530, 310, 559, 310, { color: C.amber, marker: "arrowAmber" })}
    ${lineArrow(805, 310, 834, 310, { color: C.indigo, marker: "arrowBlue" })}
    <g filter="url(#shadow)">
      <rect x="255" y="452" width="890" height="88" rx="24" fill="${C.white}" stroke="${C.teal}" stroke-opacity=".28"/>
      <rect x="255" y="452" width="12" height="88" rx="6" fill="${C.teal}"/>
      ${textBlock(285, 487, [t.formula], { size: 22, weight: 800, fill: C.ink, anchor: "start", family: "IBM Plex Mono, Noto Sans SC, ui-monospace, monospace" })}
      ${textBlock(285, 518, [t.use], { size: 15, weight: 550, fill: C.muted, anchor: "start" })}
      ${pill(1018, 477, 92, "slot 137", C.teal, C.tealSoft, { h: 34, size: 13 })}
    </g>
  `
  return frame({
    width: 1200,
    height: 575,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function sampling(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "SAMPLING",
          title: "从 Logits 到 Next Token",
          subtitle: "先约束分布，再根据温度选择确定性或随机采样路径",
          raw: "原始 Logits",
          rawBody: ["可选：先保存", "raw logprobs"],
          policy: "策略变换",
          policyBody: ["Grammar / allowlist", "penalties + temperature"],
          filter: "候选集截断",
          filterBody: ["min-p · top-k · top-p"],
          decision: "temperature",
          threshold: "< 1e-5",
          yes: "是",
          no: "否",
          greedy: "Greedy",
          greedyBody: ["argmax", "确定性选择"],
          random: "Random",
          randomBody: ["例如 Gumbel-Max", "按分布采样"],
          token: "Next-token ID",
        }
      : {
          kicker: "SAMPLING",
          title: "From Logits to the Next Token",
          subtitle:
            "Constrain the distribution first, then choose a deterministic or stochastic path",
          raw: "Raw Logits",
          rawBody: ["Optionally snapshot", "raw logprobs"],
          policy: "Policy Transform",
          policyBody: ["Grammar / allowlist", "penalties + temperature"],
          filter: "Candidate Filtering",
          filterBody: ["min-p · top-k · top-p"],
          decision: "temperature",
          threshold: "< 1e-5",
          yes: "yes",
          no: "no",
          greedy: "Greedy",
          greedyBody: ["argmax", "deterministic"],
          random: "Random",
          randomBody: ["e.g. Gumbel-Max", "sample the distribution"],
          token: "Next-token ID",
        }
  const body = `
    ${card(65, 220, 270, 135, { title: t.raw, body: t.rawBody, accent: C.blue, fill: C.blueSoft, step: "1" })}
    ${card(465, 220, 270, 135, { title: t.policy, body: t.policyBody, accent: C.indigo, fill: C.indigoSoft, step: "2" })}
    ${card(865, 220, 270, 135, { title: t.filter, body: t.filterBody, accent: C.teal, fill: C.tealSoft, step: "3" })}
    ${lineArrow(335, 287, 454, 287, { color: C.blue, marker: "arrowBlue" })}
    ${lineArrow(735, 287, 854, 287, { color: C.indigo, marker: "arrowBlue" })}
    <g filter="url(#shadow)">
      <circle cx="600" cy="450" r="76" fill="${C.white}" stroke="${C.amber}" stroke-width="2.2"/>
      ${textBlock(600, 425, [t.decision], { size: 16, weight: 700, fill: C.amber })}
      ${textBlock(600, 469, [t.threshold], { size: 25, weight: 850, fill: C.ink, family: "IBM Plex Mono, ui-monospace, monospace" })}
    </g>
    ${pathArrow("M 1000 355 C 980 410, 790 440, 684 448", { color: C.teal, marker: "arrowTeal" })}
    ${card(155, 525, 260, 115, { title: t.greedy, body: t.greedyBody, accent: C.blue, fill: C.blueSoft, align: "center" })}
    ${card(785, 525, 260, 115, { title: t.random, body: t.randomBody, accent: C.rose, fill: C.roseSoft, align: "center" })}
    ${pathArrow("M 535 490 C 450 515, 390 540, 405 560", { color: C.blue, marker: "arrowBlue" })}
    ${pathArrow("M 665 490 C 750 515, 810 540, 795 560", { color: C.rose, marker: "arrowRose" })}
    ${pill(430, 502, 58, t.yes, C.blue, C.blueSoft, { h: 28, size: 12 })}
    ${pill(712, 502, 58, t.no, C.rose, C.roseSoft, { h: 28, size: 12 })}
    <g filter="url(#shadow)">
      <rect x="490" y="592" width="220" height="60" rx="30" fill="${C.green}" />
      ${textBlock(600, 630, [t.token], { size: 19, weight: 800, fill: C.white })}
    </g>
    ${pathArrow("M 415 594 C 450 616, 470 622, 480 622", { color: C.blue, marker: "arrowBlue" })}
    ${pathArrow("M 785 594 C 750 616, 730 622, 720 622", { color: C.rose, marker: "arrowRose" })}
  `
  return frame({
    width: 1200,
    height: 690,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function preemption(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "KV PRESSURE",
          title: "请求状态与抢占循环",
          subtitle: "KV 分配失败时释放低优先级缓存，并通过重新 prefill 恢复",
          waiting: "WAITING",
          waitingBody: ["等待 token budget", "和空闲 KV block"],
          running: "RUNNING",
          runningBody: ["持续调度", "next tokens"],
          finished: "FINISHED",
          finishedBody: ["命中停止条件"],
          preempted: "PREEMPTED",
          preemptedBody: ["释放 KV blocks", "重置计算进度"],
          admit: "资源可用",
          step: "下一 step",
          stop: "停止条件",
          fail: "KV 分配失败",
          recover: "回到队列 · 重新 prefill",
        }
      : {
          kicker: "KV PRESSURE",
          title: "Request States and the Preemption Loop",
          subtitle:
            "When KV allocation fails, reclaim a low-priority cache and recover through prefill",
          waiting: "WAITING",
          waitingBody: ["Await token budget", "and free KV blocks"],
          running: "RUNNING",
          runningBody: ["Schedule the", "next tokens"],
          finished: "FINISHED",
          finishedBody: ["Stop condition reached"],
          preempted: "PREEMPTED",
          preemptedBody: ["Release KV blocks", "Reset computed progress"],
          admit: "resources available",
          step: "next step",
          stop: "stop condition",
          fail: "KV allocation fails",
          recover: "return to queue · prefill again",
        }
  const body = `
    ${card(75, 270, 260, 130, { title: t.waiting, body: t.waitingBody, accent: C.blue, fill: C.blueSoft, align: "center" })}
    ${card(470, 235, 260, 130, { title: t.running, body: t.runningBody, accent: C.teal, fill: C.tealSoft, align: "center" })}
    ${card(865, 270, 260, 130, { title: t.finished, body: t.finishedBody, accent: C.green, fill: C.greenSoft, align: "center" })}
    ${card(470, 485, 260, 130, { title: t.preempted, body: t.preemptedBody, accent: C.rose, fill: C.roseSoft, align: "center" })}
    <circle cx="44" cy="335" r="9" fill="${C.blue}"/>
    ${lineArrow(53, 335, 65, 335, { color: C.blue, marker: "arrowBlue" })}
    ${lineArrow(335, 335, 459, 300, { color: C.blue, marker: "arrowBlue" })}
    ${pill(348, 270, lang === "zh" ? 96 : 140, t.admit, C.blue, C.blueSoft, { h: 28, size: 12 })}
    ${lineArrow(730, 300, 854, 335, { color: C.green, marker: "arrowTeal" })}
    ${pill(747, 270, lang === "zh" ? 84 : 108, t.stop, C.green, C.greenSoft, { h: 28, size: 12 })}
    ${pathArrow("M 600 365 C 610 405, 610 430, 600 474", { color: C.rose, marker: "arrowRose" })}
    ${pill(625, 407, lang === "zh" ? 108 : 132, t.fail, C.rose, C.roseSoft, { h: 28, size: 12 })}
    ${pathArrow("M 470 550 C 300 560, 205 490, 205 411", { color: C.rose, marker: "arrowRose", dashed: true })}
    ${pill(198, 520, lang === "zh" ? 170 : 196, t.recover, C.rose, C.roseSoft, { h: 28, size: 12 })}
    ${pathArrow("M 535 235 C 500 180, 700 180, 665 235", { color: C.teal, marker: "arrowTeal" })}
    ${pill(555, 174, 90, t.step, C.teal, C.tealSoft, { h: 28, size: 12 })}
  `
  return frame({
    width: 1200,
    height: 665,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function gpuSystem(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "HETEROGENEOUS SYSTEM",
          title: "CPU 组织工作，GPU 并行执行",
          subtitle: "命令流与数据流通过互连进入 GPU，并在多个 SM 上展开",
          host: "HOST / CPU",
          gpu: "DEVICE / GPU",
          app: "应用程序",
          runtime: "CUDA Runtime / Driver",
          hostMem: "Host Memory",
          queue: "Command Queue",
          interconnect: "PCIe / NVLink-C2C",
          front: "GPU Front End",
          scheduler: "Global Block Scheduler",
          deviceMem: "GPU Memory + Copy Engines",
          sm: "Streaming Multiprocessors",
          warps: "resident blocks → warps",
        }
      : {
          kicker: "HETEROGENEOUS SYSTEM",
          title: "The CPU Organizes; the GPU Executes in Parallel",
          subtitle:
            "Command and data streams cross the interconnect and expand across multiple SMs",
          host: "HOST / CPU",
          gpu: "DEVICE / GPU",
          app: "Application",
          runtime: "CUDA Runtime / Driver",
          hostMem: "Host Memory",
          queue: "Command Queue",
          interconnect: "PCIe / NVLink-C2C",
          front: "GPU Front End",
          scheduler: "Global Block Scheduler",
          deviceMem: "GPU Memory + Copy Engines",
          sm: "Streaming Multiprocessors",
          warps: "resident blocks → warps",
        }
  const smCards = [0, 1, "n"]
    .map((id, index) =>
      card(885, 240 + index * 112, 245, 88, {
        title: `SM ${id}`,
        body: [t.warps],
        accent: index === 1 ? C.indigo : C.teal,
        fill: index === 1 ? C.indigoSoft : C.tealSoft,
        titleSize: 18,
        bodySize: 13,
      }),
    )
    .join("")
  const body = `
    ${section(45, 195, 365, 410, t.host, C.blue, { fill: "#f7faff" })}
    ${card(80, 275, 295, 90, { title: t.app, accent: C.blue, fill: C.blueSoft, step: "1" })}
    ${card(80, 405, 295, 90, { title: t.runtime, accent: C.indigo, fill: C.indigoSoft, step: "2" })}
    ${card(80, 535, 142, 52, { title: t.hostMem, accent: C.amber, fill: C.amberSoft, titleSize: 14 })}
    ${card(233, 535, 142, 52, { title: t.queue, accent: C.rose, fill: C.roseSoft, titleSize: 14 })}
    ${section(790, 195, 365, 410, t.gpu, C.teal, { fill: "#f6fcfa" })}
    ${card(820, 215, 205, 72, { title: t.front, accent: C.teal, fill: C.tealSoft, titleSize: 16 })}
    ${card(820, 312, 205, 72, { title: t.scheduler, accent: C.indigo, fill: C.indigoSoft, titleSize: 15 })}
    ${card(820, 505, 205, 72, { title: t.deviceMem, accent: C.amber, fill: C.amberSoft, titleSize: 14 })}
    ${smCards}
    <g filter="url(#shadow)">
      <rect x="465" y="270" width="270" height="205" rx="102" fill="${C.white}" stroke="${C.line}" stroke-width="1.5"/>
      ${textBlock(600, 340, [t.interconnect], { size: 19, weight: 800, fill: C.ink })}
      <path d="M 515 380 H 685" stroke="${C.blue}" stroke-width="9" stroke-linecap="round"/>
      <path d="M 515 412 H 685" stroke="${C.amber}" stroke-width="9" stroke-linecap="round"/>
      ${textBlock(600, 455, [lang === "zh" ? "命令流 / 数据流" : "commands / data"], { size: 14, weight: 650, fill: C.muted })}
    </g>
    ${lineArrow(375, 450, 455, 360, { color: C.blue, marker: "arrowBlue" })}
    ${lineArrow(735, 360, 810, 250, { color: C.teal, marker: "arrowTeal" })}
    ${lineArrow(222, 560, 455, 412, { color: C.amber, marker: "arrowAmber" })}
    ${lineArrow(735, 412, 810, 540, { color: C.amber, marker: "arrowAmber" })}
    ${lineArrow(1025, 348, 875, 348, { color: C.indigo, marker: "arrowBlue" })}
  `
  return frame({
    width: 1200,
    height: 650,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function kernelLaunch(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "KERNEL LAUNCH",
          title: "一次 Kernel Launch 如何抵达 SM",
          subtitle: "Host 提交通常是异步的；依赖关系或显式同步才让 CPU 等待",
          lanes: ["CPU 应用", "Runtime / Driver", "Host Queue", "GPU Front End", "SMs"],
          messages: [
            "launch kernel(grid, block, args)",
            "编码命令并入队",
            "Doorbell 通知",
            "获取并解析命令",
            "分发 Thread Blocks",
            "Blocks 完成",
          ],
          note: "CPU 可继续执行，直到遇到依赖或同步点",
        }
      : {
          kicker: "KERNEL LAUNCH",
          title: "How One Kernel Launch Reaches the SMs",
          subtitle:
            "Host submission is normally asynchronous; a dependency or explicit sync makes the CPU wait",
          lanes: ["CPU Application", "Runtime / Driver", "Host Queue", "GPU Front End", "SMs"],
          messages: [
            "launch kernel(grid, block, args)",
            "encode and enqueue work",
            "notify through a doorbell",
            "fetch and decode records",
            "distribute thread blocks",
            "blocks complete",
          ],
          note: "The CPU may continue until a dependency or synchronization point",
        }
  const xs = [125, 405, 685, 965, 1245]
  const laneColors = [C.blue, C.indigo, C.rose, C.teal, C.green]
  const lanes = xs
    .map(
      (x, i) => `<g>
        <rect x="${x - 105}" y="200" width="210" height="58" rx="18" fill="${i % 2 ? C.indigoSoft : C.blueSoft}" stroke="${laneColors[i]}" stroke-opacity=".3"/>
        ${textBlock(x, 236, [t.lanes[i]], { size: 16, weight: 800, fill: laneColors[i] })}
        <line x1="${x}" y1="270" x2="${x}" y2="610" stroke="${laneColors[i]}" stroke-opacity=".28" stroke-width="2" stroke-dasharray="7 8"/>
      </g>`,
    )
    .join("")
  const msg = (from, to, y, label, color, reverse = false) => {
    const x1 = xs[from]
    const x2 = xs[to]
    const direction = reverse
      ? "arrow"
      : color === C.teal
        ? "arrowTeal"
        : color === C.rose
          ? "arrowRose"
          : "arrowBlue"
    return `${lineArrow(x1, y, x2 + (x2 > x1 ? -8 : 8), y, { color, marker: direction, dashed: reverse })}
      ${pill(Math.min(x1, x2) + Math.abs(x2 - x1) / 2 - 110, y - 39, 220, label, color, C.white, { h: 28, size: 12 })}`
  }
  const body = `
    ${lanes}
    ${msg(0, 1, 310, t.messages[0], C.blue)}
    ${msg(1, 2, 370, t.messages[1], C.indigo)}
    ${msg(1, 3, 430, t.messages[2], C.rose)}
    ${msg(3, 2, 490, t.messages[3], C.teal, true)}
    ${msg(3, 4, 550, t.messages[4], C.teal)}
    ${msg(4, 3, 610, t.messages[5], C.green, true)}
    <g filter="url(#shadow)">
      <rect x="235" y="655" width="930" height="54" rx="27" fill="${C.ink}"/>
      ${textBlock(700, 688, [t.note], { size: 16, weight: 700, fill: C.white })}
    </g>
  `
  return frame({
    width: 1400,
    height: 750,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function programmingModel(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "PROGRAMMING MODEL",
          title: "从 Grid 到执行流水线",
          subtitle: "Block 是资源分配单位，Warp 是基本调度组",
          software: "CUDA 软件视图",
          hardware: "GPU 硬件视图",
          grid: "Grid · 一次 Kernel Launch",
          blocks: ["Thread Block 0", "Thread Block 1", "Thread Block n"],
          placed: "Block 整体放入一个 SM",
          sm: "Streaming Multiprocessor",
          warps: ["Warp 0 · 32 threads", "Warp 1 · 32 threads", "Warp n · 32 threads"],
          execution: ["CUDA Core", "Tensor Core", "Load / Store"],
        }
      : {
          kicker: "PROGRAMMING MODEL",
          title: "From a Grid to Execution Pipelines",
          subtitle: "A block is a resource-allocation unit; a warp is the basic scheduling group",
          software: "CUDA software view",
          hardware: "GPU hardware view",
          grid: "Grid · one kernel launch",
          blocks: ["Thread Block 0", "Thread Block 1", "Thread Block n"],
          placed: "A block is placed as a whole",
          sm: "Streaming Multiprocessor",
          warps: ["Warp 0 · 32 threads", "Warp 1 · 32 threads", "Warp n · 32 threads"],
          execution: ["CUDA Core", "Tensor Core", "Load / Store"],
        }
  const blocks = t.blocks
    .map((label, i) =>
      card(90, 300 + i * 92, 300, 68, {
        title: label,
        accent: C.blue,
        fill: C.blueSoft,
        titleSize: 16,
      }),
    )
    .join("")
  const warps = t.warps
    .map((label, i) =>
      card(705, 300 + i * 92, 300, 68, {
        title: label,
        accent: C.teal,
        fill: C.tealSoft,
        titleSize: 16,
      }),
    )
    .join("")
  const units = t.execution
    .map((label, i) =>
      pill(
        705 + i * 130,
        595,
        118,
        label,
        i === 1 ? C.indigo : C.amber,
        i === 1 ? C.indigoSoft : C.amberSoft,
        { h: 40, size: 12 },
      ),
    )
    .join("")
  const body = `
    ${section(45, 190, 410, 500, t.software, C.blue, { fill: "#f7faff" })}
    ${card(90, 220, 320, 62, { title: t.grid, accent: C.indigo, fill: C.indigoSoft, titleSize: 17 })}
    ${blocks}
    ${section(655, 190, 500, 500, t.hardware, C.teal, { fill: "#f6fcfa" })}
    ${card(705, 220, 400, 62, { title: t.sm, accent: C.indigo, fill: C.indigoSoft, titleSize: 17 })}
    ${warps}
    ${units}
    ${pathArrow("M 390 392 C 510 350, 560 350, 695 392", { color: C.blue, marker: "arrowBlue", width: 4 })}
    ${pill(468, 324, lang === "zh" ? 170 : 220, t.placed, C.blue, C.white, { h: 34, size: 12 })}
    ${pathArrow("M 855 576 L 855 586", { color: C.teal, marker: "arrowTeal", width: 3 })}
  `
  return frame({
    width: 1200,
    height: 735,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function divergence(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "SIMT",
          title: "Warp Divergence 如何浪费 Lane",
          subtitle: "同一 warp 的线程走不同分支时，路径分时执行，最后再汇合",
          uniform: "整齐分支",
          uniformBody: "32 个 active lanes 共同执行一条路径",
          divergent: "发散分支",
          divergentBody: "路径 A / B 分时执行，另一组 lane 被 mask",
          branch: "Warp 遇到分支",
          pathA: "执行 Path A",
          pathB: "执行 Path B",
          reconverge: "Reconverge",
          allUseful: "全部 lane 有效",
          halfMasked: "部分 lane 被 mask",
        }
      : {
          kicker: "SIMT",
          title: "How Warp Divergence Wastes Lanes",
          subtitle:
            "When threads in one warp choose different branches, paths run at different times and reconverge",
          uniform: "Uniform branch",
          uniformBody: "All 32 active lanes execute one path together",
          divergent: "Divergent branch",
          divergentBody: "Paths A and B run separately while other lanes are masked",
          branch: "Warp reaches a branch",
          pathA: "Execute Path A",
          pathB: "Execute Path B",
          reconverge: "Reconverge",
          allUseful: "all lanes useful",
          halfMasked: "some lanes masked",
        }
  const lanes = (x, y, mode, active = "all") => {
    let output = ""
    for (let i = 0; i < 32; i++) {
      const col = i % 16
      const row = Math.floor(i / 16)
      const isA = i % 2 === 0
      let fill = C.green
      let opacity = 1
      if (mode === "split") fill = isA ? C.amber : C.rose
      if (active === "a" && !isA) opacity = 0.16
      if (active === "b" && isA) opacity = 0.16
      output += `<rect x="${x + col * 18}" y="${y + row * 22}" width="14" height="16" rx="4" fill="${fill}" opacity="${opacity}"/>`
    }
    return output
  }
  const body = `
    ${card(455, 190, 290, 72, { title: t.branch, accent: C.indigo, fill: C.indigoSoft, align: "center", titleSize: 18 })}
    ${section(45, 315, 500, 310, t.uniform, C.green, { subtitle: t.allUseful, fill: C.greenSoft })}
    ${textBlock(75, 400, [t.uniformBody], { size: 15, weight: 600, fill: C.muted, anchor: "start" })}
    ${lanes(150, 435, "uniform")}
    <path d="M 120 525 H 470" stroke="${C.green}" stroke-width="12" stroke-linecap="round"/>
    ${pill(215, 548, 150, t.reconverge, C.green, C.white, { h: 34, size: 13 })}
    ${section(655, 315, 500, 310, t.divergent, C.rose, { subtitle: t.halfMasked, fill: C.roseSoft })}
    ${textBlock(685, 400, [t.divergentBody], { size: 15, weight: 600, fill: C.muted, anchor: "start" })}
    ${lanes(760, 435, "split")}
    ${pill(715, 492, 150, t.pathA, C.amber, C.amberSoft, { h: 32, size: 12 })}
    ${lanes(875, 490, "split", "a")}
    ${pill(715, 544, 150, t.pathB, C.rose, C.roseSoft, { h: 32, size: 12 })}
    ${lanes(875, 542, "split", "b")}
    ${pill(858, 592, 150, t.reconverge, C.indigo, C.white, { h: 34, size: 13 })}
    ${pathArrow("M 520 262 C 400 275, 310 290, 290 305", { color: C.green, marker: "arrowTeal" })}
    ${pathArrow("M 680 262 C 800 275, 890 290, 910 305", { color: C.rose, marker: "arrowRose" })}
  `
  return frame({
    width: 1200,
    height: 670,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function oneRecOverview(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "ONEREC",
          title: "从语义 ID 到在线反馈的完整闭环",
          subtitle: "离线表征与训练、在线生成式召回、工业过滤排序最终回流为新样本",
          offline: "离线表征与训练",
          online: "在线召回与服务",
          feedback: "反馈闭环",
          stages: [
            ["广告 / 内容表征", "RQ-VAE 量化", "多级语义 ID"],
            ["在线 / 离线样本", "序列生成模型训练", "Encoder + Decoder 发布"],
            ["用户行为与上下文", "在线特征与序列", "Beam Search 生成 SID"],
            ["SID → TID → AID", "过滤 / Quota / 粗排", "侧路精排与合并竞价"],
          ],
          impression: "曝光与反馈",
          back: "更新样本与分布",
        }
      : {
          kicker: "ONEREC",
          title: "The Full Loop from Semantic IDs to Online Feedback",
          subtitle:
            "Offline representation and training feed generative retrieval, industrial ranking, and new samples",
          offline: "Offline representation and training",
          online: "Online retrieval and serving",
          feedback: "Feedback loop",
          stages: [
            ["Ad / content representation", "RQ-VAE quantization", "Multilevel semantic ID"],
            ["Online / offline samples", "Sequence-model training", "Release Encoder + Decoder"],
            ["User behavior + context", "Online features + sequence", "Beam Search produces SIDs"],
            ["SID → TID → AID", "Filter / quota / coarse rank", "Fine rank · merge · bid"],
          ],
          impression: "Impressions and feedback",
          back: "refresh samples and distributions",
        }
  const stage = (x, y, labels, colors) =>
    labels
      .map(
        (label, i) =>
          card(x + i * 265, y, 230, 92, {
            title: label,
            accent: colors[i],
            fill: i === 0 ? C.blueSoft : i === 1 ? C.indigoSoft : C.tealSoft,
            titleSize: 15,
            align: "center",
            step: String(i + 1),
          }) +
          (i < labels.length - 1
            ? lineArrow(x + i * 265 + 230, y + 46, x + (i + 1) * 265 - 10, y + 46, {
                color: colors[i],
                marker: "arrowBlue",
              })
            : ""),
      )
      .join("")
  const body = `
    ${section(45, 190, 1110, 250, t.offline, C.indigo, { fill: "#f9f8ff" })}
    ${stage(90, 260, t.stages[0], [C.blue, C.indigo, C.teal])}
    ${stage(90, 365, t.stages[1], [C.amber, C.indigo, C.teal])}
    ${section(45, 470, 1110, 250, t.online, C.teal, { fill: "#f6fcfa" })}
    ${stage(90, 540, t.stages[2], [C.blue, C.indigo, C.teal])}
    ${stage(90, 645, t.stages[3], [C.amber, C.indigo, C.teal])}
    ${card(915, 365, 190, 92, { title: t.impression, accent: C.green, fill: C.greenSoft, align: "center", titleSize: 15 })}
    ${pathArrow("M 1010 457 C 1140 530, 1140 730, 900 750 C 600 775, 230 770, 205 468", { color: C.green, marker: "arrowTeal", dashed: true, width: 4 })}
    ${pill(488, 747, lang === "zh" ? 170 : 210, t.back, C.green, C.greenSoft, { h: 34, size: 12 })}
  `
  return frame({
    width: 1200,
    height: 810,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function retrievalPipeline(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "RETRIEVAL PIPELINE",
          title: "从入口到召回缓存的生产链路",
          subtitle: "数据接入、模型生成、创意物化和在线服务组成四个连续阶段",
          zones: ["数据接入", "模型生成", "创意物化", "在线服务"],
          steps: [
            ["入口", "Datahub"],
            ["Encoder", "Decoder"],
            ["Creative KV", "创意列表构建", "广告列表合并"],
            ["Quota / 属性过滤", "创意服务", "分片 + 粗排", "召回缓存"],
          ],
        }
      : {
          kicker: "RETRIEVAL PIPELINE",
          title: "The Production Path from Entry to Retrieval Cache",
          subtitle:
            "Data intake, model generation, creative materialization, and online serving form four stages",
          zones: ["Data intake", "Model generation", "Creative materialization", "Online serving"],
          steps: [
            ["Entry", "Datahub"],
            ["Encoder", "Decoder"],
            ["Creative KV", "Creative list build", "Ad list merge"],
            [
              "Quota / property filter",
              "Creative server",
              "Sharding + coarse score",
              "Retrieval cache",
            ],
          ],
        }
  const zoneData = [
    { x: 35, w: 240, color: C.blue, fill: "#f7faff" },
    { x: 295, w: 240, color: C.indigo, fill: "#f9f8ff" },
    { x: 555, w: 330, color: C.amber, fill: "#fffaf0" },
    { x: 905, w: 460, color: C.teal, fill: "#f6fcfa" },
  ]
  let body = ""
  let globalStep = 1
  zoneData.forEach((zone, zi) => {
    body += section(zone.x, 195, zone.w, 420, t.zones[zi], zone.color, { fill: zone.fill })
    const items = t.steps[zi]
    const spacing = (350 - items.length * 68) / Math.max(1, items.length - 1)
    items.forEach((label, i) => {
      const y = 270 + i * (68 + spacing)
      body += card(zone.x + 24, y, zone.w - 48, 64, {
        title: label,
        accent: zone.color,
        fill: C.white,
        titleSize: zi === 3 ? 14 : 15,
        step: String(globalStep++),
      })
      if (i < items.length - 1)
        body += lineArrow(zone.x + zone.w / 2, y + 64, zone.x + zone.w / 2, y + 64 + spacing - 8, {
          color: zone.color,
          marker: zone.color === C.teal ? "arrowTeal" : "arrowBlue",
        })
    })
    if (zi < zoneData.length - 1) {
      const next = zoneData[zi + 1]
      body += lineArrow(zone.x + zone.w, 405, next.x - 10, 405, {
        color: zone.color,
        marker: zone.color === C.teal ? "arrowTeal" : "arrowBlue",
        width: 4,
      })
    }
  })
  return frame({
    width: 1400,
    height: 660,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function mergePaths(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "MULTI-PATH RANKING",
          title: "主链路与 OneRec 侧路在哪里汇合",
          subtitle: "两条检索路径独立产出候选，在统一合并点去重后进入公共重排与竞价",
          mixer: "Mixer",
          main: "主链路",
          mainRetrieve: "主召回 / 粗排",
          mainRank: "主 Ranking",
          side: "OneRec 侧路",
          sideRetrieve: "GprHub 召回 / 粗排 / 缓存",
          sideRank: "Ranking Fetch + 侧路精排",
          merge: "MergeMultiAds / MergeSameAd",
          final: "公共重排与竞价",
        }
      : {
          kicker: "MULTI-PATH RANKING",
          title: "Where the Main Path and OneRec Side Path Merge",
          subtitle:
            "Two retrieval paths produce candidates independently, deduplicate once, then share reranking and bidding",
          mixer: "Mixer",
          main: "Main path",
          mainRetrieve: "Main retrieval / coarse rank",
          mainRank: "Main Ranking",
          side: "OneRec side path",
          sideRetrieve: "GprHub retrieval / coarse rank / cache",
          sideRank: "Ranking Fetch + side-path rank",
          merge: "MergeMultiAds / MergeSameAd",
          final: "Common reranking and bidding",
        }
  const body = `
    ${card(470, 190, 260, 72, { title: t.mixer, accent: C.indigo, fill: C.indigoSoft, align: "center", titleSize: 20 })}
    ${section(55, 315, 500, 230, t.main, C.blue, { fill: "#f7faff" })}
    ${card(95, 390, 190, 90, { title: t.mainRetrieve, accent: C.blue, fill: C.blueSoft, align: "center", titleSize: 15 })}
    ${card(325, 390, 190, 90, { title: t.mainRank, accent: C.indigo, fill: C.indigoSoft, align: "center", titleSize: 15 })}
    ${section(645, 315, 500, 230, t.side, C.teal, { fill: "#f6fcfa" })}
    ${card(685, 390, 190, 90, { title: t.sideRetrieve, accent: C.teal, fill: C.tealSoft, align: "center", titleSize: 14 })}
    ${card(915, 390, 190, 90, { title: t.sideRank, accent: C.amber, fill: C.amberSoft, align: "center", titleSize: 14 })}
    ${pathArrow("M 540 262 C 460 290, 260 300, 190 380", { color: C.blue, marker: "arrowBlue" })}
    ${pathArrow("M 660 262 C 740 290, 940 300, 1010 380", { color: C.teal, marker: "arrowTeal" })}
    ${lineArrow(285, 435, 315, 435, { color: C.blue, marker: "arrowBlue" })}
    ${lineArrow(875, 435, 905, 435, { color: C.teal, marker: "arrowTeal" })}
    ${card(445, 585, 310, 78, { title: t.merge, accent: C.rose, fill: C.roseSoft, align: "center", titleSize: 16 })}
    ${pathArrow("M 420 480 C 440 540, 490 560, 520 575", { color: C.blue, marker: "arrowBlue" })}
    ${pathArrow("M 1010 480 C 950 540, 780 560, 680 575", { color: C.teal, marker: "arrowTeal" })}
    ${card(445, 705, 310, 72, { title: t.final, accent: C.green, fill: C.greenSoft, align: "center", titleSize: 16 })}
    ${lineArrow(600, 663, 600, 694, { color: C.green, marker: "arrowTeal" })}
  `
  return frame({
    width: 1200,
    height: 820,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function oneRecFullOnlinePath(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "ONEREC · ONLINE TOPOLOGY",
          title: "主路与 OneRec 侧路的完整在线拓扑",
          subtitle: "线上主路与 GprHub 并发；侧路内部再选择生成式实现或旧 GPR 对照",
          request: "Mixer 请求",
          prod: "线上 / 灰度：两路并发",
          ab: "侧路内部 A/B",
          sim: "仿真隔离：可跳过主召回",
          main: "主路",
          mainCards: [
            "Retrieval Proxy",
            "过滤 · 配额 · 粗排",
            "完整 DocWash",
            "主预测 · pCTR · 出价上下文",
          ],
          side: "OneRec / GprHub 侧路",
          sideCards: [
            "生成式六节点\n或旧 GPR 对照",
            "过滤 · 配额 · 粗排",
            "异步缓存 → Ranking Fetch",
            "侧路清洗 · 预测 · pCTR",
          ],
          merge: "按 AID 合并广告 · 按创意合并与全局去重",
          rerank: "公共预测后处理 · 样式 / 模板 / 动态重排 · 拍卖输出",
          failopen: "侧路 miss / timeout / empty：主路继续",
        }
      : {
          kicker: "ONEREC · ONLINE TOPOLOGY",
          title: "The Full Online Topology of the Main and OneRec Paths",
          subtitle:
            "Production runs the main path and GprHub concurrently; the side path then selects generative or legacy GPR",
          request: "Mixer request",
          prod: "production / gray: concurrent paths",
          ab: "A/B inside the side path",
          sim: "simulation isolation: main retrieval may be skipped",
          main: "Main path",
          mainCards: [
            "Retrieval Proxy",
            "Filter · quota · coarse rank",
            "Full DocWash",
            "Main prediction · pCTR · bid context",
          ],
          side: "OneRec / GprHub side path",
          sideCards: [
            "Six generative ops\nor legacy GPR control",
            "Filter · quota · coarse rank",
            "Async cache → Ranking Fetch",
            "Side wash · prediction · pCTR",
          ],
          merge: "Merge ads by AID · merge and globally deduplicate creatives",
          rerank: "Common prediction post-process · style / template / dynamic rerank · auction",
          failopen: "side miss / timeout / empty: main path continues",
        }

  const laneCards = (x, y, labels, colors) =>
    labels
      .map((label, i) => {
        const cy = y + i * 100
        const lines = label.split("\n")
        const item = card(x, cy, 480, 76, {
          title: lines,
          accent: colors[i],
          fill: i % 2 === 0 ? C.white : colors[i] === C.teal ? C.tealSoft : C.blueSoft,
          align: "center",
          titleSize: lang === "zh" ? 16 : 15,
          step: String(i + 1),
        })
        const connector =
          i < labels.length - 1
            ? lineArrow(x + 240, cy + 76, x + 240, cy + 92, {
                color: colors[i],
                width: 2.5,
              })
            : ""
        return item + connector
      })
      .join("")

  const body = `
    ${card(555, 185, 290, 74, { title: t.request, accent: C.indigo, fill: C.indigoSoft, align: "center", titleSize: 20 })}
    ${pill(60, 282, lang === "zh" ? 215 : 270, t.prod, C.blue, C.blueSoft, { h: 32, size: 12 })}
    ${pill(555, 282, lang === "zh" ? 180 : 185, t.ab, C.teal, C.tealSoft, { h: 32, size: 12 })}
    ${pill(1015, 282, lang === "zh" ? 300 : 330, t.sim, C.amber, C.amberSoft, { h: 32, size: 12 })}
    ${section(35, 340, 620, 490, t.main, C.blue, { fill: "#f7faff" })}
    ${laneCards(105, 405, t.mainCards, [C.blue, C.blue, C.indigo, C.indigo])}
    ${section(745, 340, 620, 490, t.side, C.teal, { fill: "#f6fcfa" })}
    ${laneCards(815, 405, t.sideCards, [C.teal, C.teal, C.amber, C.teal])}
    ${pathArrow("M 625 259 C 540 300, 410 320, 345 395", { color: C.blue, width: 3 })}
    ${pathArrow("M 775 259 C 860 300, 990 320, 1055 395", { color: C.teal, width: 3 })}
    ${card(350, 885, 700, 82, { title: t.merge, accent: C.rose, fill: C.roseSoft, align: "center", titleSize: lang === "zh" ? 17 : 16 })}
    ${pathArrow("M 345 781 C 365 835, 490 855, 550 875", { color: C.blue, width: 3 })}
    ${pathArrow("M 1055 781 C 1035 835, 910 855, 850 875", { color: C.teal, width: 3 })}
    ${card(350, 1020, 700, 82, { title: t.rerank, accent: C.green, fill: C.greenSoft, align: "center", titleSize: lang === "zh" ? 16 : 15 })}
    ${lineArrow(700, 967, 700, 1009, { color: C.green, width: 3 })}
    ${pill(955, 796, 350, t.failopen, C.rose, C.white, { h: 28, size: 11 })}
    ${pathArrow("M 1295 743 C 1330 760, 1330 792, 1315 809", { color: C.rose, dashed: true, width: 2.2 })}
  `
  return frame({
    width: 1400,
    height: 1140,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function oneRecRankingMerge(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "RANKING · SOURCE WALKTHROUGH",
          title: "OneRec 在 Ranking 中怎样进入主路",
          subtitle: "两路独立完成候选绑定与预测，合并后共享重排图和最终约束",
          main: "主路候选",
          side: "GprHub 缓存候选",
          mainSteps: ["完整 DocWash", "主预测请求", "pCTR · 出价 · 生态上下文"],
          sideSteps: ["Fetch + 候选绑定", "GPR 轻量 DocWash", "侧路预测 · pCTR · 出价上下文"],
          mergeTitle: "MergeMultiAds",
          mergeBody: [
            "AID 重合：合并召回来源与策略",
            "同创意：合并模型分数；侧路独有创意追加",
            "全局创意 ID 去重：出价校正值优先，其次看新鲜度",
          ],
          post: "预测后处理与推荐上下文合并",
          rerank: "公共重排图",
          rerankBody: [["样式重排"], ["模板重排"], ["动态策略图"], ["可选 OneRec", "相关性因子"]],
          output: "统一约束 · 拍卖 · 响应",
        }
      : {
          kicker: "RANKING · SOURCE WALKTHROUGH",
          title: "How OneRec Enters the Main Path inside Ranking",
          subtitle:
            "The paths bind and score candidates independently, then share the reranking graph and final constraints",
          main: "Main-path candidates",
          side: "GprHub cached candidates",
          mainSteps: ["Full DocWash", "Main prediction request", "pCTR · bid · ecosystem context"],
          sideSteps: [
            "Fetch + candidate binding",
            "Lightweight GPR DocWash",
            "Side prediction · pCTR · bid context",
          ],
          mergeTitle: "MergeMultiAds",
          mergeBody: [
            "AID overlap: merge retrieval sources and policies",
            "Same creative: merge model scores; append side-only creatives",
            "Global creative-ID dedup: adjusted value first, then freshness",
          ],
          post: "Prediction post-process + recommendation-context merge",
          rerank: "Common reranking graph",
          rerankBody: [
            ["Style", "reranking"],
            ["Template", "reranking"],
            ["Dynamic strategy", "graph"],
            ["Optional OneRec", "relevance factor"],
          ],
          output: "Common constraints · auction · response",
        }

  const verticalSteps = (x, labels, colors) =>
    labels
      .map((label, i) => {
        const y = 300 + i * 128
        return `${card(x, y, 430, 82, { title: label, accent: colors[i], fill: i === 1 ? C.white : colors[i] === C.teal ? C.tealSoft : C.blueSoft, align: "center", titleSize: lang === "zh" ? 16 : 15, step: String(i + 1) })}${
          i < labels.length - 1
            ? lineArrow(x + 215, y + 82, x + 215, y + 116, { color: colors[i], width: 2.5 })
            : ""
        }`
      })
      .join("")

  const body = `
    ${section(45, 205, 610, 480, t.main, C.blue, { fill: "#f7faff" })}
    ${verticalSteps(135, t.mainSteps, [C.blue, C.indigo, C.blue])}
    ${section(745, 205, 610, 480, t.side, C.teal, { fill: "#f6fcfa" })}
    ${verticalSteps(835, t.sideSteps, [C.amber, C.teal, C.teal])}
    ${section(350, 735, 700, 172, t.mergeTitle, C.rose, { fill: C.roseSoft })}
    ${textBlock(700, 828, t.mergeBody, { size: lang === "zh" ? 14 : 13, weight: 600, fill: C.muted, lineHeight: 25 })}
    ${pathArrow("M 350 638 C 370 690, 470 708, 540 725", { color: C.blue, width: 3 })}
    ${pathArrow("M 1050 638 C 1025 690, 930 708, 860 725", { color: C.teal, width: 3 })}
    ${card(440, 955, 520, 70, { title: t.post, accent: C.indigo, fill: C.indigoSoft, align: "center", titleSize: lang === "zh" ? 16 : 14 })}
    ${lineArrow(700, 907, 700, 944, { color: C.indigo, width: 3 })}
    ${section(350, 1070, 700, 205, t.rerank, C.green, { fill: "#f4fbf7" })}
    ${t.rerankBody
      .map((label, i) => {
        const x = 380 + i * 162
        return card(x, 1140, 145, 82, {
          title: label,
          accent: i === 3 ? C.amber : C.green,
          fill: i === 3 ? C.amberSoft : C.white,
          align: "center",
          titleSize: lang === "zh" ? 13 : 12,
        })
      })
      .join("")}
    ${lineArrow(700, 1025, 700, 1058, { color: C.green, width: 3 })}
    ${card(500, 1320, 400, 72, { title: t.output, accent: C.green, fill: C.greenSoft, align: "center", titleSize: 17 })}
    ${lineArrow(700, 1275, 700, 1309, { color: C.green, width: 3 })}
  `
  return frame({
    width: 1400,
    height: 1440,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function debugTree(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "DIAGNOSTICS",
          title: "“最终无广告”应从哪里开始排查",
          subtitle: "沿生产链路逐级确认产物；每个缺失点都对应一组有限的责任域",
          start: "最终无广告",
          checks: [
            "Decoder 有 SID？",
            "KV 有 TID？",
            ["Creative Build", "有 AID？"],
            "公共下游保留？",
            "Fetch 命中？",
          ],
          fixes: [
            ["模型输入 / 版本", "beam"],
            ["SID 协议 / KV 版本", "分片"],
            ["TID→AID 正排", "去重"],
            ["Quota / 属性 / 库存", "粗排"],
            ["cache key / timeout", "fill_done"],
          ],
          final: ["检查 Ranking 侧路", "合并与去重"],
          yes: "是",
          no: "否",
        }
      : {
          kicker: "DIAGNOSTICS",
          title: "Where to Debug “No Final Ads”",
          subtitle:
            "Verify artifacts stage by stage; each missing product points to a bounded responsibility domain",
          start: "No final ads",
          checks: [
            ["Decoder has", "SIDs?"],
            "KV has TIDs?",
            ["Creative Build", "has AIDs?"],
            ["Common downstream", "retains ads?"],
            "Fetch hits?",
          ],
          fixes: [
            ["Model input / version", "beam"],
            ["SID protocol / KV", "version / shards"],
            ["TID→AID forward index", "dedup"],
            ["Quota / property", "inventory / coarse rank"],
            ["Cache key / timeout", "fill_done"],
          ],
          final: ["Inspect side-path ranking", "merge, and dedup"],
          yes: "yes",
          no: "no",
        }
  const colors = [C.blue, C.indigo, C.amber, C.teal, C.rose]
  let body = card(45, 250, 180, 76, {
    title: t.start,
    accent: C.rose,
    fill: C.roseSoft,
    align: "center",
    titleSize: 16,
  })
  const startX = 255
  t.checks.forEach((check, i) => {
    const x = startX + i * 220
    body += `<g filter="url(#shadow)">
      <rect x="${x}" y="220" width="190" height="116" rx="28" fill="${C.white}" stroke="${colors[i]}" stroke-width="2"/>
      <circle cx="${x + 28}" cy="248" r="17" fill="${colors[i]}"/>
      ${textBlock(x + 28, 253, [String(i + 1)], { size: 13, weight: 800, fill: C.white })}
      ${textBlock(x + 95, 287, Array.isArray(check) ? check : [check], { size: 15, weight: 750, fill: C.ink })}
    </g>`
    body += card(x, 440, 190, 118, {
      title: t.fixes[i],
      accent: colors[i],
      fill: i % 2 ? C.indigoSoft : C.blueSoft,
      align: "center",
      titleSize: 14,
    })
    body += lineArrow(x + 95, 336, x + 95, 430, {
      color: colors[i],
      marker: colors[i] === C.rose ? "arrowRose" : "arrowBlue",
      dashed: true,
    })
    body += pill(x + 108, 368, 52, t.no, colors[i], C.white, { h: 27, size: 11 })
    if (i === 0) body += lineArrow(225, 288, x - 10, 288, { color: C.rose, marker: "arrowRose" })
    if (i < t.checks.length - 1) {
      body += lineArrow(x + 190, 288, x + 210, 288, { color: colors[i], marker: "arrowBlue" })
      body += pill(x + 171, 250, 52, t.yes, colors[i], C.white, { h: 27, size: 11 })
    }
  })
  body += card(1115, 620, 240, 92, {
    title: t.final,
    accent: C.green,
    fill: C.greenSoft,
    align: "center",
    titleSize: 15,
  })
  body += pathArrow("M 1230 336 C 1275 420, 1270 530, 1235 610", {
    color: C.green,
    marker: "arrowTeal",
  })
  body += pill(1240, 410, 52, t.yes, C.green, C.white, { h: 27, size: 11 })
  return frame({
    width: 1400,
    height: 755,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function swr(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "SCHWARZ WAVEFORM RELAXATION",
          title: "在完整时间窗口上交换界面波形",
          subtitle: "每个空间子域独立求解，再用更新后的边界历史驱动下一轮迭代",
          left: "子域 1",
          right: "子域 2",
          solve: "求解完整时间窗口",
          forward: "界面波形",
          backward: "更新后的界面波形",
          time: "时间窗口",
        }
      : {
          kicker: "SCHWARZ WAVEFORM RELAXATION",
          title: "Exchange Interface Waveforms over a Full Time Window",
          subtitle:
            "Each spatial subdomain solves independently, then updated boundary histories drive the next iteration",
          left: "Subdomain 1",
          right: "Subdomain 2",
          solve: "Solve over the full time window",
          forward: "interface waveform",
          backward: "updated interface waveform",
          time: "time window",
        }
  const waveform = (x, y, color, flip = false) => {
    const d = flip
      ? `M ${x} ${y} C ${x + 35} ${y - 25}, ${x + 70} ${y + 30}, ${x + 105} ${y} S ${x + 175} ${y - 22}, ${x + 210} ${y + 4}`
      : `M ${x} ${y} C ${x + 35} ${y + 25}, ${x + 70} ${y - 30}, ${x + 105} ${y} S ${x + 175} ${y + 22}, ${x + 210} ${y - 4}`
    return `<path d="${d}" fill="none" stroke="${color}" stroke-width="5" stroke-linecap="round"/>`
  }
  const body = `
    ${section(55, 210, 420, 300, t.left, C.blue, { fill: "#f7faff" })}
    ${card(95, 275, 340, 92, { title: t.solve, accent: C.blue, fill: C.blueSoft, align: "center", titleSize: 17 })}
    ${waveform(135, 430, C.blue)}
    ${textBlock(250, 486, [t.time], { size: 13, weight: 650, fill: C.muted })}
    ${section(725, 210, 420, 300, t.right, C.teal, { fill: "#f6fcfa" })}
    ${card(765, 275, 340, 92, { title: t.solve, accent: C.teal, fill: C.tealSoft, align: "center", titleSize: 17 })}
    ${waveform(805, 430, C.teal, true)}
    ${textBlock(920, 486, [t.time], { size: 13, weight: 650, fill: C.muted })}
    ${pathArrow("M 470 285 C 560 230, 645 230, 730 285", { color: C.blue, marker: "arrowBlue", width: 4 })}
    ${pill(515, 220, 170, t.forward, C.blue, C.white, { h: 32, size: 12 })}
    ${pathArrow("M 730 455 C 645 535, 560 535, 470 455", { color: C.teal, marker: "arrowTeal", width: 4 })}
    ${pill(490, 520, 220, t.backward, C.teal, C.white, { h: 32, size: 12 })}
    <rect x="565" y="300" width="70" height="150" rx="35" fill="${C.indigoSoft}" stroke="${C.indigo}" stroke-opacity=".3"/>
    ${textBlock(600, 354, ["Γ"], { size: 38, weight: 800, fill: C.indigo })}
    ${textBlock(600, 405, [lang === "zh" ? "界面" : "interface"], { size: 13, weight: 700, fill: C.indigo })}
  `
  return frame({
    width: 1200,
    height: 610,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function methodMap(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "METHOD SELECTION",
          title: "时间并行方法选择地图",
          subtitle: "先判断线性结构，再判断空间移位求解成本或动力学耗散",
          root: "问题近似线性？",
          shifted: "复移位空间系统易解？",
          dissipative: "动力学强耗散？",
          yes: "是",
          no: "否",
          outcomes: [
            ["ParaDiag / ParaExp", "利用线性与可对角化结构"],
            ["SWR / Krylov", "使用结构感知预条件器"],
            ["Parareal / PFASST", "也可考虑 MGRiT / STMG"],
            ["特征线 / 相位校正", "谨慎使用粗时间网格"],
          ],
        }
      : {
          kicker: "METHOD SELECTION",
          title: "A Decision Map for Parallel-in-Time Methods",
          subtitle: "Test linear structure first, then spatial-shift cost or dynamical dissipation",
          root: "Approximately linear?",
          shifted: "Complex-shifted spatial systems inexpensive?",
          dissipative: "Strongly dissipative dynamics?",
          yes: "yes",
          no: "no",
          outcomes: [
            ["ParaDiag / ParaExp", "Exploit linear and diagonalizable structure"],
            ["SWR / Krylov", "Use a structure-aware preconditioner"],
            ["Parareal / PFASST", "Also consider MGRiT / STMG"],
            ["Characteristics / phase correction", "Use coarse temporal grids cautiously"],
          ],
        }
  const decision = (x, y, w, label, color) => `<g filter="url(#shadow)">
    <rect x="${x}" y="${y}" width="${w}" height="94" rx="47" fill="${C.white}" stroke="${color}" stroke-width="2.3"/>
    ${textBlock(x + w / 2, y + 58, [label], { size: 18, weight: 800, fill: C.ink })}
  </g>`
  const body = `
    ${decision(480, 190, 440, t.root, C.indigo)}
    ${decision(150, 365, 420, t.shifted, C.blue)}
    ${decision(830, 365, 420, t.dissipative, C.teal)}
    ${pathArrow("M 565 284 C 500 315, 430 325, 360 355", { color: C.blue, marker: "arrowBlue", width: 4 })}
    ${pathArrow("M 835 284 C 900 315, 970 325, 1040 355", { color: C.teal, marker: "arrowTeal", width: 4 })}
    ${pill(466, 306, 58, t.yes, C.blue, C.white, { h: 28, size: 11 })}
    ${pill(876, 306, 58, t.no, C.teal, C.white, { h: 28, size: 11 })}
    ${card(65, 560, 280, 115, { title: t.outcomes[0][0], body: [t.outcomes[0][1]], accent: C.blue, fill: C.blueSoft, align: "center", titleSize: 17, bodySize: 13 })}
    ${card(375, 560, 280, 115, { title: t.outcomes[1][0], body: [t.outcomes[1][1]], accent: C.indigo, fill: C.indigoSoft, align: "center", titleSize: 17, bodySize: 13 })}
    ${card(745, 560, 280, 115, { title: t.outcomes[2][0], body: [t.outcomes[2][1]], accent: C.teal, fill: C.tealSoft, align: "center", titleSize: 17, bodySize: 13 })}
    ${card(1055, 560, 280, 115, { title: t.outcomes[3][0], body: [t.outcomes[3][1]], accent: C.rose, fill: C.roseSoft, align: "center", titleSize: 16, bodySize: 13 })}
    ${pathArrow("M 255 459 C 230 500, 210 520, 205 550", { color: C.blue, marker: "arrowBlue" })}
    ${pathArrow("M 465 459 C 485 500, 505 520, 515 550", { color: C.indigo, marker: "arrowBlue" })}
    ${pathArrow("M 935 459 C 910 500, 890 520, 885 550", { color: C.teal, marker: "arrowTeal" })}
    ${pathArrow("M 1145 459 C 1170 500, 1190 520, 1195 550", { color: C.rose, marker: "arrowRose" })}
    ${pill(200, 492, 58, t.yes, C.blue, C.white, { h: 28, size: 11 })}
    ${pill(493, 492, 58, t.no, C.indigo, C.white, { h: 28, size: 11 })}
    ${pill(850, 492, 58, t.yes, C.teal, C.white, { h: 28, size: 11 })}
    ${pill(1142, 492, 58, t.no, C.rose, C.white, { h: 28, size: 11 })}
  `
  return frame({
    width: 1400,
    height: 730,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function sequentialTimeStepping(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "ACTA NUMERICA · FIGURE 1.1",
          title: "前向 Euler 的单向因果链",
          subtitle: "时刻 tₙ₊₁ 的状态依赖 tₙ；图中以 u₉ 为当前已知状态",
          equation: "uₙ₊₁ = uₙ + Δt f(uₙ)",
          current: "当前已知：u₉",
          pending: "尚待计算：u₁₀ 到 u₁₂",
          dependency: "先获得 u₉，才能启动 u₁₀ 的更新",
          implication: "并行方法需要重新组织跨时间段耦合，才能让多个时间区间同时工作",
        }
      : {
          kicker: "ACTA NUMERICA · FIGURE 1.1",
          title: "The One-Way Causal Chain of Forward Euler",
          subtitle: "The state at tₙ₊₁ depends on tₙ; u₉ is the latest known state in this view",
          equation: "uₙ₊₁ = uₙ + Δt f(uₙ)",
          current: "known state: u₉",
          pending: "pending: u₁₀ through u₁₂",
          dependency: "The update for u₁₀ can start only after u₉ is available",
          implication:
            "Time parallelism needs a concurrent representation of the coupling across intervals",
        }

  const xs = Array.from({ length: 13 }, (_, index) => 105 + index * 99)
  const connectors = xs
    .slice(0, -1)
    .map((x, index) => {
      const color = index < 9 ? C.teal : index === 9 ? C.rose : C.line
      return lineArrow(x + 24, 318, xs[index + 1] - 24, 318, {
        color,
        marker: color === C.teal ? "arrowTeal" : color === C.rose ? "arrowRose" : "arrow",
        dashed: index > 9,
        opacity: index > 9 ? 0.72 : 1,
      })
    })
    .join("")
  const nodes = xs
    .map((x, index) => {
      const known = index <= 9
      const active = index === 9
      const color = active ? C.rose : known ? C.teal : C.line
      const fill = active ? C.roseSoft : known ? C.tealSoft : C.slateSoft
      return `<g>
        <circle cx="${x}" cy="318" r="23" fill="${fill}" stroke="${color}" stroke-width="${active ? 3 : 2}"/>
        ${textBlock(x, 324, [`u${index}`], { size: 15, weight: 800, fill: active ? C.rose : C.ink })}
        ${textBlock(x, 370, [`t${index}`], { size: 13, weight: 650, fill: C.muted })}
      </g>`
    })
    .join("")

  const body = `
    ${card(455, 185, 490, 76, { title: t.equation, accent: C.indigo, fill: C.indigoSoft, align: "center", titleSize: 22 })}
    <path d="M 82 391 L 1320 391" fill="none" stroke="${C.line}" stroke-width="1.8" stroke-linecap="round"/>
    ${connectors}
    ${nodes}
    ${pill(780, 405, lang === "zh" ? 180 : 170, t.current, C.rose, C.roseSoft, { h: 31, size: 12 })}
    ${pill(980, 405, lang === "zh" ? 245 : 220, t.pending, C.muted, C.slateSoft, { h: 31, size: 12 })}
    ${card(85, 475, 550, 102, { title: t.dependency, accent: C.rose, fill: C.roseSoft, align: "center", titleSize: 17 })}
    ${card(765, 475, 550, 102, { title: t.implication, accent: C.blue, fill: C.blueSoft, align: "center", titleSize: 16 })}
    ${pathArrow("M 635 526 C 685 500, 715 500, 765 526", { color: C.blue, marker: "arrowBlue", width: 3 })}
  `

  return frame({
    width: 1400,
    height: 630,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function modelMemorySpectrum(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "MODEL PROBLEMS",
          title: "从时间局部性到长程记忆",
          subtitle: "扩散、边界与非线性共同决定旧信息在未来保留多久",
          left: "快速遗忘",
          right: "长期保留",
          models: [
            ["热方程", "高频先衰减", "低频与平均值"],
            ["对流扩散", "相位 + 衰减", "小黏性保留轨迹"],
            ["Burgers", "状态相关速度", "激波生成新高频"],
            ["波动方程", "双向传播与反射", "振幅和相位长存"],
          ],
          boundary: "边界会移动实际位置",
          chips: ["Dirichlet：允许出流", "Neumann：保留平均值", "周期：信号持续回返"],
        }
      : {
          kicker: "MODEL PROBLEMS",
          title: "From Temporal Locality to Long-Range Memory",
          subtitle:
            "Diffusion, boundaries, and nonlinearity determine how long past information survives",
          left: "rapid forgetting",
          right: "persistent memory",
          models: [
            ["Heat", "high frequencies decay first", "mean and slow modes remain"],
            ["Advection–diffusion", "phase + damping", "small viscosity preserves paths"],
            ["Burgers", "state-dependent velocity", "shocks regenerate fine scales"],
            ["Wave", "bidirectional propagation", "amplitude and phase persist"],
          ],
          boundary: "Boundaries shift the effective position",
          chips: [
            "Dirichlet: permits outflow",
            "Neumann: retains the mean",
            "Periodic: recirculates signals",
          ],
        }

  const colors = [C.blue, C.teal, C.amber, C.rose]
  const fills = [C.blueSoft, C.tealSoft, C.amberSoft, C.roseSoft]
  const xs = [55, 385, 715, 1045]
  const cards = t.models
    .map(
      (model, i) => `${card(xs[i], 255, 300, 185, {
        title: model[0],
        body: [model[1], model[2]],
        accent: colors[i],
        fill: fills[i],
        align: "center",
        titleSize: 21,
        bodySize: 14,
      })}
      <g opacity=".9">
        <path d="M ${xs[i] + 60} 392 C ${xs[i] + 95} ${i < 2 ? 366 : 410}, ${xs[i] + 125} ${i < 2 ? 410 : 366}, ${xs[i] + 160} 392 S ${xs[i] + 225} ${i < 2 ? 378 : 405}, ${xs[i] + 250} 392" fill="none" stroke="${colors[i]}" stroke-width="3.2" stroke-linecap="round"/>
        ${i === 0 ? `<path d="M ${xs[i] + 88} 405 L ${xs[i] + 220} 405" stroke="${colors[i]}" stroke-opacity=".25" stroke-width="2"/>` : ""}
      </g>`,
    )
    .join("")

  const body = `
    ${pill(55, 190, lang === "zh" ? 120 : 138, t.left, C.blue, C.blueSoft, { h: 30, size: 12 })}
    ${pill(1195, 190, lang === "zh" ? 150 : 150, t.right, C.rose, C.roseSoft, { h: 30, size: 12 })}
    ${pathArrow("M 185 205 C 490 170, 875 170, 1182 205", { color: C.indigo, marker: "arrowIndigo", width: 3.2 })}
    ${cards}
    ${textBlock(700, 500, [t.boundary], { size: 15, weight: 750, fill: C.indigo })}
    ${pill(130, 535, 330, t.chips[0], C.blue, C.white, { h: 42, size: 13 })}
    ${pill(535, 535, 330, t.chips[1], C.teal, C.white, { h: 42, size: 13 })}
    ${pill(940, 535, 330, t.chips[2], C.rose, C.white, { h: 42, size: 13 })}
  `
  return frame({
    width: 1400,
    height: 630,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function idcPipeline(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "IDC · PIDC · RIDC",
          title: "把残差校正排成时间流水线",
          subtitle: "流水线填满后，预测层与多个校正层在不同时间窗上同时工作",
          lanes: ["校正层 2", "校正层 1", "预测层"],
          window: "时间窗",
          fill: "启动",
          steady: "稳定并发区",
          note: "RIDC 将完整时间窗进一步缩成滑动的 M 节点窗口",
          residual: "积分残差向上游校正层传递",
        }
      : {
          kicker: "IDC · PIDC · RIDC",
          title: "Pipeline Residual Corrections across Time",
          subtitle:
            "After filling, predictor and correction levels work concurrently on different windows",
          lanes: ["Correction 2", "Correction 1", "Predictor"],
          window: "window",
          fill: "fill",
          steady: "steady concurrent region",
          note: "RIDC replaces each full window with a sliding M-node stencil",
          residual: "integral residual moves upward through correction levels",
        }
  const laneY = [245, 355, 465]
  const laneColors = [C.indigo, C.teal, C.blue]
  const laneFills = [C.indigoSoft, C.tealSoft, C.blueSoft]
  let blocks = ""
  for (let lane = 0; lane < 3; lane++) {
    const start = 2 - lane
    for (let w = 0; w < 5 - start; w++) {
      const x = 320 + (w + start) * 175
      const y = laneY[lane] - 33
      blocks += `<g filter="url(#shadow)">
        <rect x="${x}" y="${y}" width="140" height="66" rx="18" fill="${laneFills[lane]}" stroke="${laneColors[lane]}" stroke-opacity=".35"/>
        ${textBlock(x + 70, y + 40, [`${t.window} ${w}`], { size: 14, weight: 750, fill: laneColors[lane] })}
      </g>`
      if (w > 0)
        blocks += lineArrow(x - 34, laneY[lane], x - 10, laneY[lane], {
          color: laneColors[lane],
          width: 2.4,
        })
    }
  }
  let deps = ""
  for (let col = 2; col < 5; col++) {
    const x = 390 + col * 175
    deps += pathArrow(`M ${x} 432 C ${x + 18} 410, ${x + 18} 388, ${x} 378`, {
      color: C.amber,
      dashed: true,
      width: 2.1,
    })
    if (col > 2)
      deps += pathArrow(`M ${x} 322 C ${x + 18} 300, ${x + 18} 278, ${x} 268`, {
        color: C.amber,
        dashed: true,
        width: 2.1,
      })
  }
  const body = `
    ${section(55, 195, 1250, 340, lang === "zh" ? "流水线时间表" : "Pipeline schedule", C.indigo, { fill: "#fbfcff" })}
    ${laneY.map((y, i) => `${pill(78, y - 20, 180, t.lanes[i], laneColors[i], laneFills[i], { h: 40, size: 13 })}<line x1="270" y1="${y}" x2="1270" y2="${y}" stroke="${laneColors[i]}" stroke-opacity=".12" stroke-width="2"/>`).join("")}
    ${blocks}
    ${deps}
    ${pill(310, 548, 150, t.fill, C.amber, C.amberSoft, { h: 32, size: 12 })}
    ${pill(560, 548, 300, t.steady, C.green, C.greenSoft, { h: 32, size: 12 })}
    ${pill(920, 548, 390, t.residual, C.amber, C.amberSoft, { h: 32, size: 11 })}
    ${card(190, 600, 1020, 72, { title: t.note, accent: C.teal, fill: C.tealSoft, align: "center", titleSize: 15 })}
  `
  return frame({
    width: 1400,
    height: 710,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function paraExpDecomposition(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "PARAEXP",
          title: "局部受迫响应与齐次尾部分开计算",
          subtitle: "线性叠加把各时间段的零初值解重构为完整演化",
          intervals: ["区间 1", "区间 2", "区间 3", "区间 4"],
          local: "零初值受迫解",
          tail: "矩阵指数传播",
          parallel: "第一步：四个 vₙ 并行",
          reconstruction: "第二步：叠加到目标区间",
          formula: "u(t) = vⱼ(t) + Σₙ<ⱼ wₙ(t)",
        }
      : {
          kicker: "PARAEXP",
          title: "Separate Local Forced Responses from Homogeneous Tails",
          subtitle:
            "Linearity reconstructs the full evolution from zero-initial-value interval solves",
          intervals: ["Interval 1", "Interval 2", "Interval 3", "Interval 4"],
          local: "zero-initial forced solve",
          tail: "matrix-exponential tail",
          parallel: "Step 1: all vₙ solves run concurrently",
          reconstruction: "Step 2: superpose contributions on the target interval",
          formula: "u(t) = vⱼ(t) + Σₙ<ⱼ wₙ(t)",
        }
  const xs = [75, 375, 675, 975]
  const colors = [C.blue, C.teal, C.indigo, C.rose]
  const fills = [C.blueSoft, C.tealSoft, C.indigoSoft, C.roseSoft]
  const cards = xs
    .map((x, i) =>
      card(x, 240, 250, 145, {
        title: `v${i + 1}(t)`,
        body: [t.intervals[i], t.local],
        accent: colors[i],
        fill: fills[i],
        align: "center",
        titleSize: 22,
        bodySize: 13,
      }),
    )
    .join("")
  const tails = [0, 1, 2]
    .map((i) => {
      const x1 = xs[i] + 230
      const x2 = xs[3] + 125
      const y = 425 + i * 46
      return `${pathArrow(`M ${x1} 385 C ${x1 + 70} ${y}, ${x2 - 120} ${y}, ${x2 - 45} ${y}`, { color: colors[i], width: 2.5 })}${pill(x1 + 35, y - 17, 150, `w${i + 1}(t)`, colors[i], C.white, { h: 30, size: 12 })}`
    })
    .join("")
  const body = `
    ${pill(75, 190, lang === "zh" ? 250 : 300, t.parallel, C.blue, C.blueSoft, { h: 34, size: 12 })}
    ${cards}
    ${tails}
    ${pill(930, 405, 275, t.tail, C.rose, C.roseSoft, { h: 34, size: 12 })}
    ${card(360, 560, 680, 90, { title: t.formula, body: [t.reconstruction], accent: C.green, fill: C.greenSoft, align: "center", titleSize: 22, bodySize: 13 })}
    ${pathArrow("M 1110 515 C 1030 545, 980 555, 1040 590", { color: C.green, width: 2.8 })}
  `
  return frame({
    width: 1400,
    height: 690,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function paraDiagStages(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "PARADIAG",
          title: "时间变换把全时间系统拆成独立空间求解",
          subtitle: "变换与逆变换是全局步骤，中间的复移位空间系统可以完全并发",
          input: "全时间系统",
          transform: "时间方向变换",
          solve: "独立移位空间系统",
          inverse: "逆时间变换",
          output: "全部时间状态",
          parallel: "并行区域",
          notes: ["FFT / V⁻¹", "Dⱼ ⊗ Iₓ − A", "FFT⁻¹ / V"],
        }
      : {
          kicker: "PARADIAG",
          title: "A Time Transform Splits the All-at-Once System",
          subtitle:
            "The transforms are global; all complex-shifted spatial systems in the middle are concurrent",
          input: "All-at-once system",
          transform: "Transform in time",
          solve: "Shifted spatial systems",
          inverse: "Inverse time transform",
          output: "states",
          parallel: "concurrent region",
          notes: ["FFT / V⁻¹", "Dⱼ ⊗ Iₓ − A", "FFT⁻¹ / V"],
        }
  let systems = ""
  for (let i = 0; i < 4; i++) {
    systems += `<g><rect x="610" y="${230 + i * 70}" width="300" height="52" rx="16" fill="${i % 2 ? C.tealSoft : C.indigoSoft}" stroke="${i % 2 ? C.teal : C.indigo}" stroke-opacity=".35"/>${textBlock(760, 263 + i * 70, [`shift ${i + 1}: (d${i + 1}I − A)x = b`], { size: 14, weight: 700, fill: C.ink })}</g>`
  }
  const body = `
    ${card(55, 285, 210, 130, { title: t.input, body: ["K U = b"], accent: C.blue, fill: C.blueSoft, align: "center", titleSize: 18, bodySize: 20 })}
    ${card(330, 285, 210, 130, { title: t.transform, body: [t.notes[0]], accent: C.amber, fill: C.amberSoft, align: "center", titleSize: 17, bodySize: 15 })}
    ${section(575, 190, 370, 390, t.solve, C.indigo, { subtitle: t.parallel, fill: "#fbfaff" })}
    ${systems}
    ${card(980, 285, 210, 130, { title: t.inverse, body: [t.notes[2]], accent: C.amber, fill: C.amberSoft, align: "center", titleSize: 17, bodySize: 15 })}
    ${card(1250, 285, 100, 130, { title: "U", body: [t.output], accent: C.green, fill: C.greenSoft, align: "center", titleSize: 24, bodySize: 12 })}
    ${lineArrow(265, 350, 320, 350, { color: C.blue, width: 3 })}
    ${lineArrow(540, 350, 565, 350, { color: C.amber, width: 3 })}
    ${lineArrow(945, 350, 970, 350, { color: C.indigo, width: 3 })}
    ${lineArrow(1190, 350, 1240, 350, { color: C.amber, width: 3 })}
    ${pill(655, 595, 210, t.notes[1], C.indigo, C.indigoSoft, { h: 36, size: 13 })}
  `
  return frame({
    width: 1400,
    height: 665,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function parabolicMultilevelMap(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "PARABOLIC PINT",
          title: "四种方法如何利用粗层传递慢模态",
          subtitle: "扩散先削弱局部高频，粗时间或粗时空表示再承担长程校正",
          cards: [
            ["Parareal", "并行细传播", "顺序粗预测"],
            ["PFASST", "配置节点 SDC", "FAS 粗 sweep"],
            ["MGRiT", "F / C 点松弛", "递归时间粗层"],
            ["STMG", "时间块 Jacobi", "空间 + 时间粗化"],
          ],
          top: "局部细节：并行处理",
          bottom: "慢变分量：粗层长程传递",
          warning: "黏性下降后，相位误差进入粗层，四种方法都会逐步退化",
        }
      : {
          kicker: "PARABOLIC PINT",
          title: "How Four Methods Carry Slow Modes on Coarse Levels",
          subtitle:
            "Diffusion first suppresses local high frequencies; coarse time or space–time levels then communicate globally",
          cards: [
            ["Parareal", "concurrent fine propagation", "sequential coarse prediction"],
            ["PFASST", "SDC on collocation nodes", "coarse FAS sweep"],
            ["MGRiT", "F/C relaxation", "recursive temporal levels"],
            ["STMG", "time-block Jacobi", "space + time coarsening"],
          ],
          top: "local detail: concurrent work",
          bottom: "slow components: long-range coarse correction",
          warning:
            "As viscosity falls, phase error enters the coarse level and all four methods degrade",
        }
  const xs = [55, 385, 715, 1045]
  const colors = [C.blue, C.teal, C.indigo, C.green]
  const fills = [C.blueSoft, C.tealSoft, C.indigoSoft, C.greenSoft]
  const cards = t.cards
    .map(
      (m, i) => `${section(xs[i], 220, 300, 315, m[0], colors[i], { fill: C.white })}
    ${card(xs[i] + 25, 290, 250, 88, { title: m[1], accent: colors[i], fill: fills[i], align: "center", titleSize: 15 })}
    ${pathArrow(`M ${xs[i] + 150} 378 C ${xs[i] + 120} 405, ${xs[i] + 120} 425, ${xs[i] + 150} 445`, { color: colors[i], width: 2.4 })}
    ${card(xs[i] + 25, 445, 250, 66, { title: m[2], accent: C.amber, fill: C.amberSoft, align: "center", titleSize: 14 })}`,
    )
    .join("")
  const body = `
    ${pill(55, 175, 300, t.top, C.blue, C.blueSoft, { h: 32, size: 12 })}
    ${pill(1000, 175, 345, t.bottom, C.amber, C.amberSoft, { h: 32, size: 12 })}
    ${cards}
    ${card(210, 580, 980, 72, { title: t.warning, accent: C.rose, fill: C.roseSoft, align: "center", titleSize: 15 })}
  `
  return frame({
    width: 1400,
    height: 690,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function longSequenceDesignSpace(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "长序列推荐",
          title: "五条路线把成本约束放在不同位置",
          subtitle: "从可复用的用户状态，逐步走向候选相关的事件级交互",
          routes: [
            ["压缩", "L → M 个兴趣", "预算键：M", "长期状态可复用"],
            ["层次化", "局部块 → 全局 token", "预算键：块长 B", "保留时间结构"],
            ["线性模型", "序列 → 递归状态", "预算键：状态维度 d", "长度近似线性扩展"],
            ["检索", "L → Top-Kᵣ 行为", "预算键：Kᵣ", "候选相关性优先"],
            ["稀疏交互", "每点保留 Kₛ 条边", "预算键：Kₛ / W", "保留更多原始细节"],
          ],
          left: "可复用的用户状态",
          right: "候选相关的精细交互",
          dense: "稠密端到端：信息路径最完整，计算与显存成本最高",
        }
      : {
          kicker: "LONG-SEQUENCE REC",
          title: "Five Routes Place the Cost Constraint at Different Points",
          subtitle:
            "The design spectrum moves from reusable user state to candidate-specific event interaction",
          routes: [
            ["Compression", "L → M interests", "budget key: M", "reusable long-term state"],
            [
              "Hierarchy",
              "local blocks → global tokens",
              "budget key: block B",
              "preserves temporal structure",
            ],
            [
              "Linear model",
              "sequence → recurrent state",
              "budget key: state width d",
              "near-linear length scaling",
            ],
            ["Retrieval", "L → Top-Kᵣ events", "budget key: Kᵣ", "prioritizes target relevance"],
            [
              "Sparse interaction",
              "Kₛ edges per event",
              "budget key: Kₛ / W",
              "retains more raw detail",
            ],
          ],
          left: "reusable user state",
          right: "candidate-specific interaction",
          dense:
            "Dense end to end: the most complete paths, with the highest compute and memory cost",
        }
  const xs = [40, 310, 580, 850, 1120]
  const colors = [C.teal, C.green, C.indigo, C.blue, C.amber]
  const fills = [C.tealSoft, C.greenSoft, C.indigoSoft, C.blueSoft, C.amberSoft]
  const cards = t.routes
    .map(
      (route, i) => `${section(xs[i], 215, 240, 300, route[0], colors[i], { fill: C.white })}
      ${pill(xs[i] + 18, 287, 204, route[1], colors[i], fills[i], { h: 38, size: 13 })}
      ${textBlock(xs[i] + 120, 376, [route[2]], { size: 15, weight: 750, fill: colors[i] })}
      ${textBlock(xs[i] + 120, 428, [route[3]], { size: 14, weight: 550, fill: C.muted })}
      <circle cx="${xs[i] + 120}" cy="478" r="7" fill="${colors[i]}" opacity=".8"/>`,
    )
    .join("")
  const body = `
    ${cards}
    ${pathArrow("M 78 575 C 390 612, 1010 612, 1322 575", { color: C.teal, width: 3 })}
    ${textBlock(58, 646, [t.left], { size: 14, weight: 700, fill: C.teal, anchor: "start" })}
    ${textBlock(1342, 646, [t.right], { size: 14, weight: 700, fill: C.amber, anchor: "end" })}
    ${pill(295, 670, 810, t.dense, C.rose, C.roseSoft, { h: 36, size: 13 })}
  `
  return frame({
    width: 1400,
    height: 750,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function longSequenceHybridMemory(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "混合架构",
          title: "长期压缩记忆与候选精细路径并行工作",
          subtitle: "两条路径使用不同容量与刷新频率，在融合层汇合",
          history: "完整行为历史",
          historyBody: ["长度 L · 长期偏好"],
          compressor: "兴趣压缩器",
          compressorBody: ["学习查询 · 异步更新"],
          interests: "M 个兴趣 token",
          interestsBody: ["跨候选复用 · 固定容量"],
          candidate: "当前候选 q",
          candidateBody: ["ID · 上下文 · 多模态"],
          sparse: "语义检索 + 局部窗口",
          sparseBody: ["Top-Kᵣ 行为 / 每点 Kₛ 条边"],
          events: "候选相关行为子图",
          eventsBody: ["最新事件 · 保留原始细节"],
          fusion: "目标感知融合",
          fusionBody: ["长期状态 + 精细行为"],
          output: "排序表示 / 分数",
          memory: "长期记忆路径",
          detail: "候选精细路径",
          fallback: "语义索引异常时：近期窗口 + 缓存兴趣",
        }
      : {
          kicker: "HYBRID ARCHITECTURE",
          title: "Compressed Long-Term Memory and Candidate Detail Run in Parallel",
          subtitle:
            "The paths use different capacities and refresh rates, then meet at target-aware fusion",
          history: "Full behavior history",
          historyBody: ["length L · long-term signals"],
          compressor: "Interest compressor",
          compressorBody: ["learned queries · async refresh"],
          interests: "M interest tokens",
          interestsBody: ["reusable · fixed capacity"],
          candidate: "Current candidate q",
          candidateBody: ["ID · context · multimodal"],
          sparse: "Semantic sparse selection",
          sparseBody: ["Top-Kᵣ events / Kₛ edges"],
          events: "Candidate-relevant subgraph",
          eventsBody: ["latest events · raw detail"],
          fusion: "Target-aware fusion",
          fusionBody: ["long-term state + event detail"],
          output: "Ranking state / score",
          memory: "long-term memory path",
          detail: "candidate-detail path",
          fallback: "Index failure fallback: recent window + cached interests",
        }
  const body = `
    ${section(35, 195, 1030, 205, t.memory, C.teal, { fill: "#fbfffe" })}
    ${card(65, 260, 230, 105, { title: t.history, body: t.historyBody, accent: C.teal, fill: C.tealSoft, titleSize: 16, bodySize: 12 })}
    ${card(380, 260, 250, 105, { title: t.compressor, body: t.compressorBody, accent: C.green, fill: C.greenSoft, titleSize: 16, bodySize: 12 })}
    ${card(720, 260, 260, 105, { title: t.interests, body: t.interestsBody, accent: C.indigo, fill: C.indigoSoft, titleSize: 16, bodySize: 12 })}
    ${lineArrow(295, 312, 368, 312, { color: C.teal, width: 2.8 })}
    ${lineArrow(630, 312, 708, 312, { color: C.green, width: 2.8 })}

    ${section(35, 435, 1030, 205, t.detail, C.blue, { fill: "#fbfdff" })}
    ${card(65, 500, 230, 105, { title: t.candidate, body: t.candidateBody, accent: C.blue, fill: C.blueSoft, titleSize: 16, bodySize: 12 })}
    ${card(380, 500, 250, 105, { title: t.sparse, body: t.sparseBody, accent: C.amber, fill: C.amberSoft, titleSize: 15, bodySize: 12 })}
    ${card(720, 500, 260, 105, { title: t.events, body: t.eventsBody, accent: C.rose, fill: C.roseSoft, titleSize: 15, bodySize: 12 })}
    ${lineArrow(295, 552, 368, 552, { color: C.blue, width: 2.8 })}
    ${lineArrow(630, 552, 708, 552, { color: C.amber, width: 2.8 })}

    ${card(1120, 310, 225, 190, { title: t.fusion, body: t.fusionBody, accent: C.indigo, fill: C.white, align: "center", titleSize: 17, bodySize: 12 })}
    ${pathArrow("M 980 312 C 1050 312, 1050 362, 1108 385", { color: C.indigo, width: 3 })}
    ${pathArrow("M 980 552 C 1050 552, 1050 470, 1108 435", { color: C.rose, width: 3 })}
    ${pathArrow("M 295 552 C 790 720, 1110 670, 1228 512", { color: C.blue, width: 2.1, dashed: true, opacity: 0.75 })}
    ${pill(1088, 555, 265, t.output, C.green, C.greenSoft, { h: 38, size: 14 })}
    ${lineArrow(1232, 500, 1232, 543, { color: C.green, width: 2.5 })}
    ${pill(480, 677, 440, t.fallback, C.amber, C.amberSoft, { h: 36, size: 13 })}
  `
  return frame({
    width: 1400,
    height: 750,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function longSequenceDecisionFramework(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "技术选型",
          title: "先识别主要信号，再检查工程门槛",
          subtitle: "每个场景给出起点；最终方案仍需在等成本条件下做消融",
          question: "哪一种信号决定收益？",
          cases: [
            ["目标相关性", "广告 / 强关联", "检索", "Top-Kᵣ"],
            ["极长 + 多模态", "内容流 / 兴趣漂移", "稀疏 + 压缩", "Kₛ + M"],
            ["稳定长期偏好", "电商 / 工具", "压缩", "M"],
            ["局部时间连续", "会话 / 周期行为", "层次化", "B + W"],
            ["显存或长度约束", "研究基线 / 离线", "线性模型", "d"],
          ],
          gate: "工程门槛",
          checks: [
            ["索引与表示", "版本、刷新、陈旧率"],
            ["运行时", "真稀疏 kernel、低精度、批处理"],
            ["可靠性", "近期窗口 / 缓存兴趣回退"],
            ["证据", "相同硬件、预算、输入分布"],
          ],
        }
      : {
          kicker: "METHOD SELECTION",
          title: "Identify the Dominant Signal, Then Check Engineering Gates",
          subtitle:
            "Each setting suggests a starting point; the final choice still needs equal-cost ablation",
          question: "Which signal drives value?",
          cases: [
            ["Target relevance", "ads / tight relation", "Retrieval", "Top-Kᵣ"],
            ["Extreme + multimodal", "feeds / rapid drift", "Sparse + compression", "Kₛ + M"],
            ["Stable long-term taste", "commerce / utility", "Compression", "M"],
            ["Local temporal continuity", "sessions / periodicity", "Hierarchy", "B + W"],
            ["Memory or length limit", "research baseline / offline", "Linear model", "d"],
          ],
          gate: "engineering gates",
          checks: [
            ["Index + representation", "versions, refresh, staleness"],
            ["Runtime", "true sparse kernels, low precision, batching"],
            ["Reliability", "recent-window / cached-interest fallback"],
            ["Evidence", "same hardware, budget, and input distribution"],
          ],
        }
  const xs = [30, 305, 580, 855, 1130]
  const colors = [C.blue, C.amber, C.teal, C.green, C.indigo]
  const fills = [C.blueSoft, C.amberSoft, C.tealSoft, C.greenSoft, C.indigoSoft]
  const cases = t.cases
    .map(
      (entry, i) => `${section(xs[i], 250, 240, 250, entry[0], colors[i], { fill: C.white })}
      ${textBlock(xs[i] + 120, 332, [entry[1]], { size: 14, weight: 550, fill: C.muted })}
      ${pill(xs[i] + 20, 372, 200, entry[2], colors[i], fills[i], { h: 44, size: 15 })}
      ${textBlock(xs[i] + 120, 462, [entry[3]], { size: 17, weight: 800, fill: colors[i] })}`,
    )
    .join("")
  const checkXs = [65, 395, 725, 1055]
  const checks = t.checks
    .map(
      (entry, i) =>
        `${card(checkXs[i], 585, 280, 90, { title: entry[0], body: [entry[1]], accent: [C.teal, C.blue, C.rose, C.green][i], fill: C.white, titleSize: 15, bodySize: 12 })}`,
    )
    .join("")
  const body = `
    ${pill(505, 185, 390, t.question, C.indigo, C.indigoSoft, { h: 40, size: 16 })}
    ${cases}
    ${textBlock(54, 548, [t.gate.toUpperCase()], { size: 13, weight: 800, fill: C.muted, anchor: "start" })}
    <path d="M 55 560 L 1345 560" stroke="${C.line}" stroke-width="1.5" stroke-dasharray="5 8"/>
    ${checks}
  `
  return frame({
    width: 1400,
    height: 720,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

const diagrams = [
  ["ml-inference", "serving-loop", servingLoop],
  ["ml-inference", "paged-kv-cache", pagedKv],
  ["ml-inference", "sampling-pipeline", sampling],
  ["ml-inference", "preemption-states", preemption],
  ["gpu", "heterogeneous-system", gpuSystem],
  ["gpu", "kernel-launch", kernelLaunch],
  ["gpu", "programming-model", programmingModel],
  ["gpu", "warp-divergence", divergence],
  ["onerec", "overview", oneRecOverview],
  ["onerec", "retrieval-pipeline", retrievalPipeline],
  ["onerec", "merge-paths", mergePaths],
  ["onerec", "full-online-path", oneRecFullOnlinePath],
  ["onerec", "ranking-merge-detail", oneRecRankingMerge],
  ["onerec", "debugging-tree", debugTree],
  ["pint", "schwarz-waveform-relaxation", swr],
  ["pint", "method-selection", methodMap],
  ["pint", "sequential-time-stepping", sequentialTimeStepping],
  ["pint", "model-memory-spectrum", modelMemorySpectrum],
  ["pint", "idc-pipeline", idcPipeline],
  ["pint", "paraexp-decomposition", paraExpDecomposition],
  ["pint", "paradiag-three-stage", paraDiagStages],
  ["pint", "parabolic-multilevel-map", parabolicMultilevelMap],
  ["long-sequence-recommendation", "design-space", longSequenceDesignSpace],
  ["long-sequence-recommendation", "hybrid-memory", longSequenceHybridMemory],
  ["long-sequence-recommendation", "decision-framework", longSequenceDecisionFramework],
]

let generated = 0
for (const [group, name, builder] of diagrams) {
  for (const lang of ["zh", "en"]) {
    const outputDir = path.join(assetRoot, group, lang)
    fs.mkdirSync(outputDir, { recursive: true })
    const svg = builder(lang).replace(/[ \t]+$/gm, "")
    fs.writeFileSync(path.join(outputDir, `${name}.svg`), svg, "utf8")
    generated++
  }
}

console.log(`Generated ${generated} SVG diagrams in ${path.relative(repoRoot, assetRoot)}`)
