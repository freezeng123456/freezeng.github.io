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

const publications = JSON.parse(
  fs.readFileSync(path.join(scriptDir, "data", "tao-zhou-publications.json"), "utf8"),
)

function researchMap(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "RESEARCH MAP",
          title: "七个专题共享同一条技术主线",
          subtitle: "先把不确定性写成可逼近的对象，再从数据反推参数，最后让时间演化可算",
          groups: [
            "把随机输入变成可逼近的对象",
            "从有限数据反推未知量",
            "让时间演化在有限资源下可算",
          ],
          themes: [
            [
              "随机逼近与配点设计",
              "23 篇",
              "把 Galerkin、配点与最小二乘",
              "统一成加权采样问题",
              "核心量：Christoffel 函数",
            ],
            [
              "谱方法与降阶表示",
              "6 篇",
              "无界区域基函数、低秩流形",
              "与分数幂矩阵的快速算法",
              "核心量：解的衰减与秩",
            ],
            [
              "贝叶斯反问题与数据同化",
              "10 篇",
              "代理模型与采样器交替细化",
              "误差指标决定何时补算真解",
              "核心量：后验误差指标",
            ],
            [
              "科学机器学习",
              "22 篇",
              "残差损失、密度流与算子学习",
              "采样点分布是主要自由度",
              "核心量：失效概率与频谱",
            ],
            [
              "FBSDE 与随机最优控制",
              "22 篇",
              "用概率表示替代网格离散",
              "多步插值提高时间精度",
              "核心量：条件期望",
            ],
            [
              "相场模型与变步长离散",
              "12 篇",
              "变步长 BDF 与 IMEX 的",
              "离散能量与二次型正定性",
              "核心量：步长比阈值",
            ],
            [
              "时间并行算法",
              "13 篇",
              "全时间系统的对角化与预条件",
              "把串行递推变成并发求解",
              "核心量：谱与条件数",
            ],
          ],
          moves: [
            "把串行或高维结构写成一个整体算子",
            "为该算子设计带权、带自适应的近似逆",
            "用误差指标决定在哪里补算真实信息",
          ],
        }
      : {
          kicker: "RESEARCH MAP",
          title: "Seven Topics on One Technical Spine",
          subtitle:
            "First make randomness approximable, then invert data for unknowns, then make time evolution computable",
          groups: [
            "turn random input into approximable objects",
            "recover unknowns from limited data",
            "make time evolution affordable",
          ],
          themes: [
            [
              "Stochastic approximation",
              "23 papers",
              "Galerkin, collocation and least squares",
              "become one weighted sampling problem",
              "central object: Christoffel function",
            ],
            [
              "Spectral and reduced order",
              "6 papers",
              "bases on unbounded domains, low-rank",
              "manifolds, fractional matrix powers",
              "central object: decay and rank",
            ],
            [
              "Bayesian inverse problems",
              "10 papers",
              "surrogate and sampler refine together",
              "an indicator decides when to solve",
              "central object: posterior indicator",
            ],
            [
              "Scientific machine learning",
              "22 papers",
              "residual losses, density flows, operators",
              "the sample distribution is the knob",
              "central object: failure probability",
            ],
            [
              "FBSDEs and control",
              "22 papers",
              "probabilistic representation instead",
              "of grids; multistep interpolation",
              "central object: conditional expectation",
            ],
            [
              "Phase field and time steps",
              "12 papers",
              "discrete energy and quadratic-form",
              "positivity for variable-step BDF/IMEX",
              "central object: step-ratio threshold",
            ],
            [
              "Parallel in time",
              "13 papers",
              "diagonalise and precondition the",
              "all-at-once space-time operator",
              "central object: spectrum, conditioning",
            ],
          ],
          moves: [
            "write the sequential or high-dimensional structure as one operator",
            "design a weighted, adaptive approximate inverse for it",
            "let an error indicator decide where true information is added",
          ],
        }

  const groupX = [30, 496, 962]
  const groupW = [420, 420, 408]
  const groupColor = [C.blue, C.indigo, C.teal]
  const themeColor = [C.blue, C.green, C.indigo, C.rose, C.teal, C.amber, C.blue]
  const themeFill = [
    C.blueSoft,
    C.greenSoft,
    C.indigoSoft,
    C.roseSoft,
    C.tealSoft,
    C.amberSoft,
    C.blueSoft,
  ]
  const layout = [
    [0, 1],
    [2, 3],
    [4, 5, 6],
  ]

  let panels = ""
  layout.forEach((indices, group) => {
    panels += section(groupX[group], 196, groupW[group], 466, t.groups[group], groupColor[group], {
      fill: C.white,
    })
    const count = indices.length
    const height = count === 3 ? 122 : 190
    const gap = count === 3 ? 132 : 200
    indices.forEach((themeIndex, slot) => {
      const theme = t.themes[themeIndex]
      const y = 262 + slot * gap
      panels += card(groupX[group] + 22, y, groupW[group] - 44, height, {
        title: theme[0],
        body: count === 3 ? [theme[2], theme[3]] : [theme[2], theme[3], theme[4]],
        accent: themeColor[themeIndex],
        fill: themeFill[themeIndex],
        titleSize: 18,
        bodySize: 13,
      })
      panels += pill(
        groupX[group] + groupW[group] - 128,
        y + 14,
        84,
        theme[1],
        themeColor[themeIndex],
        C.white,
        { h: 26, size: 12 },
      )
    })
  })

  const body = `
    ${panels}
    ${lineArrow(456, 429, 490, 429, { color: C.indigo, width: 2.8 })}
    ${lineArrow(922, 429, 956, 429, { color: C.teal, width: 2.8 })}
    ${textBlock(700, 700, [lang === "zh" ? "三个反复出现的动作" : "Three recurring moves"], {
      size: 16,
      weight: 800,
      fill: C.ink,
    })}
    ${pill(30, 726, 420, t.moves[0], C.blue, C.white, { h: 44, size: 12 })}
    ${pill(496, 726, 420, t.moves[1], C.indigo, C.white, { h: 44, size: 12 })}
    ${pill(962, 726, 408, t.moves[2], C.teal, C.white, { h: 44, size: 12 })}
  `
  return frame({
    width: 1400,
    height: 800,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function researchTimeline(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "TIMELINE",
          title: "2010 至 2026 年的专题分布",
          subtitle: "圆点面积表示该年在该专题上的论文数量；数据取自主页发表列表",
          lanes: {
            "stochastic-approximation": "随机逼近与配点设计",
            "spectral-reduced-order": "谱方法与降阶表示",
            "bayesian-inference": "贝叶斯反问题",
            "scientific-machine-learning": "科学机器学习",
            "fbsde-control": "FBSDE 与随机控制",
            "phase-field": "相场与变步长离散",
            "parallel-in-time": "时间并行算法",
          },
          legend: "1 篇 / 2 篇 / 4 篇及以上",
          note: "论文按发表年份计入；投稿与预印本按主页标注年份计入",
        }
      : {
          kicker: "TIMELINE",
          title: "Topic Distribution from 2010 to 2026",
          subtitle: "Disc area encodes papers per topic per year, taken from the homepage list",
          lanes: {
            "stochastic-approximation": "Stochastic approximation",
            "spectral-reduced-order": "Spectral and reduced order",
            "bayesian-inference": "Bayesian inverse problems",
            "scientific-machine-learning": "Scientific machine learning",
            "fbsde-control": "FBSDEs and control",
            "phase-field": "Phase field, time steps",
            "parallel-in-time": "Parallel in time",
          },
          legend: "1 / 2 / 4 or more papers",
          note: "Counted by publication year; submissions use the year given on the homepage",
        }

  const laneOrder = [
    "stochastic-approximation",
    "spectral-reduced-order",
    "bayesian-inference",
    "scientific-machine-learning",
    "fbsde-control",
    "phase-field",
    "parallel-in-time",
  ]
  const laneColor = [C.blue, C.green, C.indigo, C.rose, C.teal, C.amber, C.blue]
  const years = Array.from({ length: 17 }, (_, index) => 2010 + index)
  const left = 330
  const step = 62
  const x = (year) => left + (year - 2010) * step

  const counts = new Map()
  for (const paper of publications.papers) {
    const key = `${paper.theme}:${paper.year}`
    counts.set(key, (counts.get(key) ?? 0) + 1)
  }

  let grid = years
    .map(
      (year) =>
        `<path d="M ${x(year)} 216 L ${x(year)} 640" stroke="${C.line}" stroke-opacity=".35" stroke-width="1"/>
         ${textBlock(x(year), 668, [`${year}`], { size: 12, weight: 650, fill: C.muted })}`,
    )
    .join("")

  let lanes = ""
  laneOrder.forEach((theme, index) => {
    const y = 250 + index * 56
    lanes += `<rect x="30" y="${y - 24}" width="1340" height="48" rx="16" fill="${index % 2 ? C.slateSoft : C.white}" fill-opacity=".7"/>
      ${textBlock(300, y + 5, [t.lanes[theme]], { size: 14, weight: 700, fill: laneColor[index], anchor: "end" })}`
    for (const year of years) {
      const count = counts.get(`${theme}:${year}`) ?? 0
      if (!count) continue
      const radius = 7 + Math.sqrt(count) * 5
      lanes += `<circle cx="${x(year)}" cy="${y}" r="${radius}" fill="${laneColor[index]}" fill-opacity=".2" stroke="${laneColor[index]}" stroke-width="1.8"/>
        ${textBlock(x(year), y + 5, [`${count}`], { size: 13, weight: 800, fill: laneColor[index] })}`
    }
  })

  const body = `
    ${grid}
    ${lanes}
    <circle cx="80" cy="716" r="12" fill="${C.blue}" fill-opacity=".2" stroke="${C.blue}" stroke-width="1.8"/>
    <circle cx="124" cy="716" r="14" fill="${C.blue}" fill-opacity=".2" stroke="${C.blue}" stroke-width="1.8"/>
    <circle cx="174" cy="716" r="17" fill="${C.blue}" fill-opacity=".2" stroke="${C.blue}" stroke-width="1.8"/>
    ${textBlock(204, 721, [t.legend], { size: 13, weight: 650, fill: C.muted, anchor: "start" })}
    ${pill(760, 698, lang === "zh" ? 500 : 560, t.note, C.muted, C.white, { h: 36, size: 12 })}
  `
  return frame({
    width: 1400,
    height: 770,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function samplingDesign(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "WEIGHTED LEAST SQUARES",
          title: "配点设计的统一流程",
          subtitle: "同一条流水线容纳随机采样、诱导采样与贪心选点；差别只在“从哪个密度取点”",
          stages: [
            ["选定多项式空间", "指标集 Λ，维数 N", "正交基 {φₙ}", "目标测度 ω"],
            ["构造采样密度", "μ 与 ω 的比值", "决定权函数 w", "Christoffel 函数 λ_N"],
            ["取 M 个样本点", "随机 / 诱导 / 贪心", "形成设计矩阵 A", "权重矩阵 W"],
            ["解加权最小二乘", "min ‖W(Ac − f)‖₂", "或稀疏化 ℓ₁ 问题", "得到系数 c"],
          ],
          gramian: "关键量：加权 Gram 矩阵",
          gramianBody: "样本足够多时它接近单位矩阵，最小二乘才稳定",
          budget: "样本预算随采样密度改变",
          budgets: [
            ["均匀 / 目标测度采样", "M ≳ N² 量级", C.rose],
            ["Christoffel 加权采样", "M ≳ N log N", C.amber],
            ["贪心选点（近似 Fekete）", "M 接近 N", C.green],
          ],
          feedback: "同一份候选点集可反复复用，只更换权与选点准则",
        }
      : {
          kicker: "WEIGHTED LEAST SQUARES",
          title: "One Pipeline for Collocation Design",
          subtitle:
            "Random, induced and greedy designs share this pipeline; only the sampling density changes",
          stages: [
            [
              "Fix the polynomial space",
              "index set Λ, dimension N",
              "basis {φₙ}",
              "target measure ω",
            ],
            [
              "Build a sampling density",
              "ratio of μ against ω",
              "fixes the weight w",
              "Christoffel function λ_N",
            ],
            [
              "Draw M design points",
              "random / induced / greedy",
              "assemble the matrix A",
              "weight matrix W",
            ],
            [
              "Weighted least squares",
              "min ‖W(Ac − f)‖₂",
              "or a sparse ℓ₁ problem",
              "return coefficients c",
            ],
          ],
          gramian: "Central object: the weighted Gram matrix",
          gramianBody: "Least squares is stable once enough samples make it close to the identity",
          budget: "The sample budget follows the sampling density",
          budgets: [
            ["uniform / target-measure draws", "M of order N²", C.rose],
            ["Christoffel-weighted draws", "M of order N log N", C.amber],
            ["greedy points (approximate Fekete)", "M close to N", C.green],
          ],
          feedback: "One candidate pool is reused; only weights and the selection rule change",
        }

  const xs = [40, 385, 730, 1075]
  const colors = [C.blue, C.indigo, C.teal, C.green]
  const fills = [C.blueSoft, C.indigoSoft, C.tealSoft, C.greenSoft]
  const stages = t.stages
    .map(
      (stage, i) => `${card(xs[i], 200, 285, 196, {
        title: stage[0],
        body: [stage[1], stage[2], stage[3]],
        accent: colors[i],
        fill: fills[i],
        step: `${i + 1}`,
        titleSize: 17,
        bodySize: 13,
      })}
      ${i < 3 ? lineArrow(xs[i] + 292, 298, xs[i + 1] - 8, 298, { color: colors[i], width: 3 }) : ""}`,
    )
    .join("")

  const budgets = t.budgets
    .map(
      (row, i) => `<g>
        ${pill(742, 500 + i * 58, 380, row[0], row[2], C.white, { h: 44, size: 13 })}
        ${pill(1140, 500 + i * 58, 225, row[1], row[2], C.white, { h: 44, size: 13 })}
      </g>`,
    )
    .join("")

  const body = `
    ${stages}
    ${card(40, 470, 640, 200, {
      title: t.gramian,
      body: ["G = Aᵀ W² A ≈ I_N", t.gramianBody],
      accent: C.amber,
      fill: C.amberSoft,
      align: "center",
      titleSize: 18,
      bodySize: 15,
    })}
    ${textBlock(1055, 470, [t.budget], { size: 16, weight: 800, fill: C.ink })}
    ${budgets}
    ${pill(lang === "zh" ? 420 : 400, 420, lang === "zh" ? 560 : 600, t.feedback, C.muted, C.white, { h: 34, size: 13 })}
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

function sparseRecovery(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "SPARSE RECOVERY",
          title: "把梯度信息与采样密度一起写进恢复问题",
          subtitle: "同一个稀疏系数向量同时解释函数值与导数值，测量行数因此成倍增加",
          valueRow: "函数值测量",
          gradRow: "梯度测量",
          valueBody: "每个样本 1 行",
          valueFormula: "φ(y⁽ᵐ⁾)ᵀ c ≈ f(y⁽ᵐ⁾)",
          gradBody: "同一样本再加 d 行",
          gradFormula: "∂φ(y⁽ᵐ⁾)ᵀ c ≈ ∂f(y⁽ᵐ⁾)",
          stack: "堆叠测量矩阵",
          stackBody: ["行数由 M 变为 M(1+d)", "未知量仍是 N 个系数"],
          precond: "预条件与采样密度",
          precondBody: ["左乘对角权矩阵", "使各行范数接近一致"],
          rip: "有界正交系下的样本复杂度",
          ripBody: [
            "行归一化后落入有界正交系框架",
            "s 稀疏恢复所需行数约为 s 乘对数因子",
            "行范数越均衡，常数越小",
          ],
          program: "ℓ₁ 或加权 ℓ₁ 问题",
          programBody: ["min ‖c‖₁", "s.t. ‖W(Ac − b)‖₂ ≤ δ"],
          gain: "每个昂贵样本的信息量",
          gainBody: [
            "一次正问题求解同时给出函数值",
            "与全部偏导数，行数放大 (1+d) 倍",
            "未知量个数保持不变",
          ],
        }
      : {
          kicker: "SPARSE RECOVERY",
          title: "Fold Gradients and Sampling Density into the Recovery Problem",
          subtitle:
            "One sparse coefficient vector must explain values and derivatives, multiplying the measurement rows",
          valueRow: "value measurements",
          gradRow: "gradient measurements",
          valueBody: "one row per sample",
          valueFormula: "φ(y⁽ᵐ⁾)ᵀ c ≈ f(y⁽ᵐ⁾)",
          gradBody: "d more rows per sample",
          gradFormula: "∂φ(y⁽ᵐ⁾)ᵀ c ≈ ∂f(y⁽ᵐ⁾)",
          stack: "Stacked measurement matrix",
          stackBody: ["rows grow from M to M(1+d)", "while N coefficients stay unknown"],
          precond: "Preconditioning and density",
          precondBody: ["a diagonal weight matrix", "equalises the row norms"],
          rip: "Sample complexity",
          ripBody: [
            "row scaling puts the matrix in the",
            "bounded orthonormal system setting,",
            "where s-sparse recovery needs O(s · log)",
          ],
          program: "ℓ₁ or weighted ℓ₁ program",
          programBody: ["min ‖c‖₁", "s.t. ‖W(Ac − b)‖₂ ≤ δ"],
          gain: "Information per expensive sample",
          gainBody: [
            "one forward solve returns the value",
            "and every partial derivative, so rows",
            "grow by (1+d) with no new unknowns",
          ],
        }

  const rows = (x, y, w, color, count) => {
    let out = ""
    for (let i = 0; i < count; i++) {
      out += `<rect x="${x}" y="${y + i * 17}" width="${w}" height="11" rx="5" fill="${color}" fill-opacity="${0.22 + 0.1 * (i % 3)}"/>`
    }
    return out
  }

  const body = `
    ${section(40, 196, 400, 252, t.valueRow, C.blue, { fill: "#f7faff", subtitle: t.valueBody })}
    ${rows(72, 262, 336, C.blue, 4)}
    ${pill(80, 358, 320, t.valueFormula, C.blue, C.white, { h: 38, size: 14 })}
    ${textBlock(240, 424, [lang === "zh" ? "共 M 行、N 列" : "M rows, N columns"], { size: 14, weight: 700, fill: C.blue })}
    ${section(40, 476, 400, 252, t.gradRow, C.amber, { fill: "#fffaf0", subtitle: t.gradBody })}
    ${rows(72, 542, 336, C.amber, 6)}
    ${pill(80, 650, 320, t.gradFormula, C.amber, C.white, { h: 38, size: 14 })}
    ${textBlock(240, 706, [lang === "zh" ? "共 Md 行、N 列" : "Md rows, N columns"], { size: 14, weight: 700, fill: C.amber })}
    ${pathArrow("M 448 300 C 480 300, 480 372, 496 372", { color: C.blue, width: 3 })}
    ${pathArrow("M 448 600 C 480 600, 480 460, 496 460", { color: C.amber, width: 3 })}
    ${card(504, 340, 390, 156, {
      title: t.stack,
      body: t.stackBody,
      accent: C.indigo,
      fill: C.indigoSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${card(504, 196, 390, 128, {
      title: t.precond,
      body: t.precondBody,
      accent: C.teal,
      fill: C.tealSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${lineArrow(699, 330, 699, 334, { color: C.teal, width: 3 })}
    ${card(504, 552, 390, 148, {
      title: t.program,
      body: t.programBody,
      accent: C.green,
      fill: C.greenSoft,
      align: "center",
      titleSize: 17,
      bodySize: 15,
    })}
    ${lineArrow(699, 502, 699, 546, { color: C.indigo, width: 3 })}
    ${section(950, 196, 410, 300, t.rip, C.rose, { fill: "#fff8f9" })}
    ${textBlock(1155, 360, t.ripBody, { size: 14, weight: 550, fill: C.muted, lineHeight: 30 })}
    ${lineArrow(902, 400, 942, 400, { color: C.rose, width: 3 })}
    ${section(950, 528, 410, 200, t.gain, C.blue, { fill: "#f8fbff" })}
    ${textBlock(1155, 640, t.gainBody, { size: 14, weight: 550, fill: C.muted, lineHeight: 30 })}
    ${lineArrow(902, 626, 942, 626, { color: C.green, width: 3 })}
  `
  return frame({
    width: 1400,
    height: 770,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function bayesianSurrogateLoop(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "ADAPTIVE INFERENCE",
          title: "代理模型与采样器交替细化",
          subtitle: "只有当误差指标判定当前代理不可信时，才调用一次昂贵的正问题求解",
          nodes: [
            ["先验与观测", "高斯先验", "加噪声模型", "确定后验密度"],
            ["当前代理", "多项式混沌、神经网络", "或算子网络", "给出快速前向映射"],
            ["采样或反演", "MCMC、RTO、SVGD", "或集合 Kalman", "都在代理上推进"],
            ["误差指标", "残差、方差", "或校正量大小", "判断代理是否可信"],
            ["补算真解", "在被选中的点上", "调用真前向模型", "并扩充训练集"],
          ],
          cheap: "廉价循环：代理内部反复采样",
          expensive: "昂贵步骤：真前向模型",
          criterion: "触发准则",
          criterionBody: "指标超过阈值即细化；否则继续在代理上采样",
          guard: "去偏保护",
          guardBody: "接受率修正或后验加权，防止代理误差直接进入后验",
        }
      : {
          kicker: "ADAPTIVE INFERENCE",
          title: "Surrogate and Sampler Refine Together",
          subtitle:
            "An expensive forward solve happens only when the error indicator declares the surrogate untrustworthy",
          nodes: [
            ["Prior and data", "Gaussian prior", "and a noise model", "fix the posterior"],
            [
              "Current surrogate",
              "polynomial chaos,",
              "network or operator net",
              "fast forward map",
            ],
            ["Sample or invert", "MCMC, RTO, SVGD", "or ensemble Kalman", "run on the surrogate"],
            ["Error indicator", "residual, variance", "or correction size", "judges the surrogate"],
            ["True solve", "the exact forward model", "at selected points", "extends the data set"],
          ],
          cheap: "cheap loop: repeated surrogate evaluations",
          expensive: "expensive step: the true forward model",
          criterion: "Refinement trigger",
          criterionBody: "refine when the indicator exceeds a threshold, otherwise keep sampling",
          guard: "Debiasing guard",
          guardBody:
            "acceptance corrections or posterior reweighting keep surrogate error out of the posterior",
        }

  const nx = [31, 303, 575, 847, 1119]
  const colors = [C.blue, C.indigo, C.teal, C.amber, C.rose]
  const fills = [C.blueSoft, C.indigoSoft, C.tealSoft, C.amberSoft, C.roseSoft]
  const nodes = t.nodes
    .map(
      (node, i) => `${card(nx[i], 250, 250, 150, {
        title: node[0],
        body: [node[1], node[2], node[3]],
        accent: colors[i],
        fill: fills[i],
        step: `${i + 1}`,
        titleSize: 16,
        bodySize: 12,
      })}
      ${i < 4 ? lineArrow(nx[i] + 257, 325, nx[i + 1] - 8, 325, { color: colors[i], width: 3 }) : ""}`,
    )
    .join("")

  const body = `
    ${nodes}
    ${pathArrow("M 1244 406 C 1244 476, 430 476, 430 408", { color: C.rose, width: 3 })}
    ${pathArrow("M 640 244 C 640 208, 372 208, 372 242", { color: C.teal, width: 2.6, dashed: true })}
    ${pill(lang === "zh" ? 590 : 550, 458, lang === "zh" ? 400 : 470, t.expensive, C.rose, C.white, { h: 34, size: 12 })}
    ${pill(lang === "zh" ? 400 : 370, 194, lang === "zh" ? 300 : 360, t.cheap, C.teal, C.white, { h: 30, size: 12 })}
    ${card(31, 550, 660, 160, {
      title: t.criterion,
      body: [t.criterionBody],
      accent: C.amber,
      fill: C.amberSoft,
      align: "center",
      titleSize: 18,
      bodySize: 14,
    })}
    ${card(709, 550, 660, 160, {
      title: t.guard,
      body: [t.guardBody],
      accent: C.green,
      fill: C.greenSoft,
      align: "center",
      titleSize: 18,
      bodySize: 14,
    })}
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

function failureInformedSampling(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "FAILURE-INFORMED PINN",
          title: "用可靠性分析决定下一批配点",
          subtitle: "把“残差偏大”定义成失效事件，再用失效概率控制采样与停机",
          field: "残差场",
          fieldNote: "训练后的 PINN 在部分区域残差明显偏大",
          event: "失效事件",
          eventBody: "定义 {y : ε(y) > 阈值}，其中 ε 为局部残差",
          prob: "失效概率",
          probBody: "P_F = ∫ 1{ε(y)>阈值} dρ(y)，用采样估计",
          sample: "在失效域内取新点",
          sampleBody: "重要采样或子集模拟给出失效域内的样本",
          retrain: "扩充训练集并继续训练",
          retrainBody: "旧点保留或按重采样策略替换",
          stop: "停机准则",
          stopBody: "失效概率降到容差以下即停止，不再依赖固定迭代次数",
          inverse: "反问题变体",
          inverseBody: "观测残差与 PDE 残差共同定义失效事件",
        }
      : {
          kicker: "FAILURE-INFORMED PINN",
          title: "Let Reliability Analysis Choose the Next Collocation Batch",
          subtitle:
            "Declare a large residual to be a failure event, then use its probability to drive sampling and stopping",
          field: "residual field",
          fieldNote: "after training, some regions keep a visibly larger residual",
          event: "Failure event",
          eventBody: "define {y : ε(y) > tolerance} with ε the local residual",
          prob: "Failure probability",
          probBody: "P_F = ∫ 1{ε(y)>tol} dρ(y), estimated by sampling",
          sample: "Draw points inside the failure region",
          sampleBody: "importance sampling or subset simulation supplies them",
          retrain: "Extend the training set and continue",
          retrainBody: "old points are kept or replaced by a resampling rule",
          stop: "Stopping rule",
          stopBody:
            "stop when the failure probability drops below tolerance, not after a fixed count",
          inverse: "Inverse-problem variant",
          inverseBody: "data misfit and PDE residual jointly define the failure event",
        }

  let heat = ""
  for (let i = 0; i < 8; i++) {
    for (let j = 0; j < 6; j++) {
      const hot = (i >= 5 && j >= 3) || (i <= 1 && j <= 1)
      heat += `<rect x="${72 + i * 44}" y="${248 + j * 30}" width="40" height="26" rx="7" fill="${hot ? C.rose : C.blue}" fill-opacity="${hot ? 0.32 + 0.08 * ((i + j) % 3) : 0.08}"/>`
    }
  }

  const body = `
    ${section(40, 196, 400, 300, t.field, C.blue, { fill: C.white, subtitle: "ε(y)" })}
    ${heat}
    ${textBlock(240, 466, [t.fieldNote], { size: 13, weight: 600, fill: C.muted })}
    ${card(470, 200, 420, 132, {
      title: t.event,
      body: [t.eventBody],
      accent: C.rose,
      fill: C.roseSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${card(470, 356, 420, 132, {
      title: t.prob,
      body: [t.probBody],
      accent: C.amber,
      fill: C.amberSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${card(920, 200, 440, 132, {
      title: t.sample,
      body: [t.sampleBody],
      accent: C.indigo,
      fill: C.indigoSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${card(920, 356, 440, 132, {
      title: t.retrain,
      body: [t.retrainBody],
      accent: C.teal,
      fill: C.tealSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${lineArrow(447, 266, 462, 266, { color: C.rose, width: 3 })}
    ${lineArrow(680, 338, 680, 350, { color: C.rose, width: 3 })}
    ${lineArrow(897, 422, 912, 422, { color: C.amber, width: 3 })}
    ${lineArrow(1140, 338, 1140, 350, { color: C.indigo, width: 3 })}
    ${pathArrow("M 1000 494 C 1000 552, 240 552, 240 510", { color: C.teal, width: 3, dashed: true })}
    ${card(40, 580, 640, 150, {
      title: t.stop,
      body: [t.stopBody],
      accent: C.green,
      fill: C.greenSoft,
      align: "center",
      titleSize: 17,
      bodySize: 14,
    })}
    ${card(720, 580, 640, 150, {
      title: t.inverse,
      body: [t.inverseBody],
      accent: C.blue,
      fill: C.blueSoft,
      align: "center",
      titleSize: 17,
      bodySize: 14,
    })}
  `
  return frame({
    width: 1400,
    height: 770,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function densityFlowSolvers(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "FLOW-BASED PDE SOLVERS",
          title: "用可逆映射把密度方程写成参数化问题",
          subtitle: "解不再是网格上的数值，而是一个把参考密度推前到目标密度的可逆变换",
          refTitle: "参考密度",
          refBody: ["标准高斯或均匀分布", "可直接采样"],
          mapTitle: "可逆映射 f_θ",
          mapBody: ["由若干耦合层复合", "Jacobian 为三角结构", "行列式可显式计算"],
          pushTitle: "推前密度",
          pushBody: ["p(x) = q(f_θ(x)) · |det ∇f_θ(x)|", "q 是参考密度，无需归一化常数"],
          lossTitle: "损失来自方程本身",
          lossBody: ["把 p_θ 代入 Fokker-Planck 残差", "在采样点上取平方平均"],
          timeTitle: "时间方向",
          timeBody: ["把 t 作为条件变量输入", "一次训练覆盖整个时间区间"],
          boundedTitle: "有界支撑的处理",
          boundedBody: ["逻辑映射把无界坐标压回区间", "避免密度泄漏到定义域之外"],
          gain: "两个自动满足的性质",
          gainList: ["非负性由推前公式保证", "归一化由变量替换保证"],
        }
      : {
          kicker: "FLOW-BASED PDE SOLVERS",
          title: "Write a Density Equation as a Parametrised Transport Map",
          subtitle:
            "The unknown becomes an invertible map that pushes a reference density onto the target",
          refTitle: "Reference density",
          refBody: ["standard Gaussian or uniform", "can be sampled directly"],
          mapTitle: "Invertible map f_θ",
          mapBody: [
            "a composition of coupling layers",
            "triangular Jacobian",
            "explicit determinant",
          ],
          pushTitle: "Push-forward density",
          pushBody: [
            "p(x) = q(f_θ(x)) · |det ∇f_θ(x)|",
            "q is the reference; no normaliser needed",
          ],
          lossTitle: "The loss comes from the equation",
          lossBody: [
            "insert p_θ into the Fokker-Planck residual",
            "average its square over samples",
          ],
          timeTitle: "The time direction",
          timeBody: ["feed t as a conditioning input", "one training run covers the interval"],
          boundedTitle: "Handling bounded supports",
          boundedBody: [
            "a logistic map folds coordinates back",
            "so no mass leaks outside the domain",
          ],
          gain: "Two properties hold by construction",
          gainList: [
            "positivity follows from the push-forward",
            "normalisation follows from the change of variables",
          ],
        }

  let layers = ""
  for (let i = 0; i < 4; i++) {
    const x = 520 + i * 92
    layers += `<rect x="${x}" y="252" width="72" height="164" rx="16" fill="${i % 2 ? C.indigoSoft : C.tealSoft}" stroke="${i % 2 ? C.indigo : C.teal}" stroke-opacity=".4"/>
      ${textBlock(x + 36, 342, [`L${i + 1}`], { size: 17, weight: 800, fill: i % 2 ? C.indigo : C.teal })}`
  }

  const body = `
    ${card(40, 250, 240, 166, {
      title: t.refTitle,
      body: t.refBody,
      accent: C.blue,
      fill: C.blueSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${section(500, 196, 400, 268, t.mapTitle, C.indigo, { fill: C.white })}
    ${layers}
    ${textBlock(700, 442, [t.mapBody[1]], { size: 13, weight: 650, fill: C.muted })}
    ${lineArrow(288, 333, 492, 333, { color: C.blue, width: 3 })}
    ${card(940, 250, 420, 166, {
      title: t.pushTitle,
      body: t.pushBody,
      accent: C.teal,
      fill: C.tealSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${lineArrow(908, 333, 932, 333, { color: C.indigo, width: 3 })}
    ${card(40, 496, 420, 180, {
      title: t.lossTitle,
      body: t.lossBody,
      accent: C.rose,
      fill: C.roseSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${card(490, 496, 420, 180, {
      title: t.timeTitle,
      body: t.timeBody,
      accent: C.amber,
      fill: C.amberSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${card(940, 496, 420, 180, {
      title: t.boundedTitle,
      body: t.boundedBody,
      accent: C.green,
      fill: C.greenSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${pathArrow("M 1150 424 C 1150 460, 260 460, 260 490", { color: C.rose, width: 2.6, dashed: true })}
    ${textBlock(700, 716, [t.gain], { size: 15, weight: 800, fill: C.ink })}
    ${pill(300, 736, lang === "zh" ? 340 : 400, t.gainList[0], C.teal, C.white, { h: 34, size: 12 })}
    ${pill(lang === "zh" ? 680 : 740, 736, lang === "zh" ? 400 : 440, t.gainList[1], C.blue, C.white, { h: 34, size: 12 })}
  `
  return frame({
    width: 1400,
    height: 800,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function operatorLearningUq(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "OPERATOR LEARNING + UQ",
          title: "算子学习中的三处不确定性",
          subtitle: "输入函数的采样方式、算子近似本身与外推区域各自贡献不同的误差来源",
          inputTitle: "输入函数样本",
          inputBody: ["每个样本是一条离散化的函数", "传感器位置可以不固定"],
          encoderTitle: "编码器",
          encoderBody: ["置换不变的集合编码", "或固定网格分支网络"],
          decoderTitle: "解码器",
          decoderBody: ["在查询点上求值", "与坐标嵌入结合"],
          predTitle: "预测输出",
          predBody: ["每个查询点给出一个值", "以及一个置信区间"],
          uqTitle: "不确定性来源",
          rows: [
            ["输入采样", "传感器数量与位置", C.blue],
            ["算子近似", "有限参数与有限数据", C.indigo],
            ["外推", "查询落在训练分布之外", C.rose],
          ],
          methods: "三条给出置信度的路线",
          methodList: [
            ["信息瓶颈", "压缩输入表示，用变分界给出噪声尺度"],
            ["隐变量 + 高斯过程", "隐变量给均值，过程给相关结构"],
            ["集合或后验采样", "多组参数给出预测分布"],
          ],
        }
      : {
          kicker: "OPERATOR LEARNING + UQ",
          title: "Three Places Uncertainty Enters Operator Learning",
          subtitle:
            "How input functions are sampled, how the operator is approximated, and where queries fall",
          inputTitle: "Input function samples",
          inputBody: [
            "each sample is a discretised function",
            "sensor locations need not be fixed",
          ],
          encoderTitle: "Encoder",
          encoderBody: ["permutation-invariant set encoding", "or a fixed-grid branch network"],
          decoderTitle: "Decoder",
          decoderBody: ["evaluates at query points", "combined with coordinate embeddings"],
          predTitle: "Prediction",
          predBody: ["a value at every query point", "together with a confidence band"],
          uqTitle: "Sources of uncertainty",
          rows: [
            ["Input sampling", "sensor count and placement", C.blue],
            ["Operator approximation", "finite parameters, finite data", C.indigo],
            ["Extrapolation", "queries outside the training law", C.rose],
          ],
          methods: "Three routes to a calibrated confidence",
          methodList: [
            ["Information bottleneck", "compress the input, bound the noise variationally"],
            ["Latent variable and GP", "the latent gives a mean, the process a covariance"],
            ["Ensembles or posteriors", "several parameter sets give a predictive spread"],
          ],
        }

  let curves = ""
  for (let i = 0; i < 4; i++) {
    const y = 262 + i * 38
    curves += `<path d="M 72 ${y} C 118 ${y - 16 + i * 4}, 168 ${y + 18 - i * 5}, 214 ${y - 4} S 300 ${y + 12}, 340 ${y - 8}" fill="none" stroke="${C.blue}" stroke-width="2.6" stroke-opacity="${0.35 + i * 0.15}" stroke-linecap="round"/>`
  }

  const rows = t.rows
    .map(
      (row, i) => `${pill(720, 500 + i * 60, 300, row[0], row[2], C.white, { h: 46, size: 14 })}
      ${pill(1036, 500 + i * 60, 330, row[1], row[2], C.white, { h: 46, size: 13 })}`,
    )
    .join("")

  const methods = t.methodList
    .map(
      (
        method,
        i,
      ) => `${pill(40, 500 + i * 60, lang === "zh" ? 190 : 250, method[0], [C.teal, C.green, C.amber][i], C.white, { h: 46, size: 13 })}
      ${pill(lang === "zh" ? 246 : 306, 500 + i * 60, lang === "zh" ? 400 : 340, method[1], [C.teal, C.green, C.amber][i], C.white, { h: 46, size: 12 })}`,
    )
    .join("")

  const body = `
    ${section(40, 196, 340, 232, t.inputTitle, C.blue, { fill: C.white })}
    ${curves}
    ${textBlock(210, 408, [t.inputBody[1]], { size: 12, weight: 600, fill: C.muted })}
    ${card(420, 214, 290, 196, {
      title: t.encoderTitle,
      body: t.encoderBody,
      accent: C.indigo,
      fill: C.indigoSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${card(750, 214, 290, 196, {
      title: t.decoderTitle,
      body: t.decoderBody,
      accent: C.teal,
      fill: C.tealSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${card(1080, 214, 280, 196, {
      title: t.predTitle,
      body: t.predBody,
      accent: C.rose,
      fill: C.roseSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${lineArrow(388, 312, 412, 312, { color: C.blue, width: 3 })}
    ${lineArrow(718, 312, 742, 312, { color: C.indigo, width: 3 })}
    ${lineArrow(1048, 312, 1072, 312, { color: C.teal, width: 3 })}
    ${textBlock(340, 468, [t.methods], { size: 15, weight: 800, fill: C.ink })}
    ${textBlock(1050, 468, [t.uqTitle], { size: 15, weight: 800, fill: C.ink })}
    ${methods}
    ${rows}
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

function fbsdeMultistep(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "FBSDE SCHEMES",
          title: "多步格式如何提高倒向方程的时间精度",
          subtitle: "倒推一步时使用多个已知的未来层，用插值多项式代替被积函数",
          forward: "正向：模拟状态轨道",
          backward: "倒向：解条件期望",
          known: "已知层",
          target: "待求层",
          interpTitle: "多步插值",
          interp: "用 k 个已知未来层做 Lagrange 插值",
          zproc: "Z 过程由同一插值给出",
          quadTitle: "条件期望的求积",
          quad: "对 Brown 增量用 Gauss-Hermite 求积",
          space: "空间方向用插值或稀疏网格表示",
          orderTitle: "阶数与代价",
          order: "k 个未来层给出约 k 阶时间精度",
          cost: "代价是启动值与更严的稳定性条件",
        }
      : {
          kicker: "FBSDE SCHEMES",
          title: "How Multistep Schemes Raise the Backward Accuracy",
          subtitle:
            "Each backward step uses several known future levels, replacing the integrand by an interpolant",
          forward: "forward: simulate the state paths",
          backward: "backward: solve conditional expectations",
          known: "known levels",
          target: "level being solved",
          interpTitle: "Multistep interpolation",
          interp: "Lagrange interpolation over k future levels",
          zproc: "the Z process uses the same interpolant",
          quadTitle: "Quadrature for expectations",
          quad: "Gauss-Hermite rule in the Brownian increment",
          space: "space uses interpolation or a sparse grid",
          orderTitle: "Order and cost",
          order: "k future levels give roughly order k in time",
          cost: "the price is starting values and stability",
        }

  const xs = Array.from({ length: 8 }, (_, i) => 120 + i * 152)
  let axis = `<path d="M 80 300 L 1340 300" stroke="${C.line}" stroke-width="1.8" stroke-linecap="round"/>
    <path d="M 80 470 L 1340 470" stroke="${C.line}" stroke-width="1.8" stroke-linecap="round"/>`
  let nodes = ""
  xs.forEach((x, i) => {
    nodes += `<circle cx="${x}" cy="300" r="18" fill="${C.tealSoft}" stroke="${C.teal}" stroke-width="2"/>
      ${textBlock(x, 306, [`X${i}`], { size: 13, weight: 750, fill: C.teal })}`
    const isTarget = i === 3
    const isKnown = i > 3
    const color = isTarget ? C.rose : isKnown ? C.indigo : C.line
    const fill = isTarget ? C.roseSoft : isKnown ? C.indigoSoft : C.slateSoft
    nodes += `<circle cx="${x}" cy="470" r="18" fill="${fill}" stroke="${color}" stroke-width="${isTarget ? 3 : 2}"/>
      ${textBlock(x, 476, [`Y${i}`], { size: 13, weight: 750, fill: color === C.line ? C.muted : color })}`
    if (i < 7) {
      nodes += lineArrow(x + 24, 300, xs[i + 1] - 24, 300, { color: C.teal, width: 2.4 })
    }
    if (i > 3) {
      nodes += `<path d="M ${x} 448 C ${x - 20} 400, ${xs[3] + 30} 400, ${xs[3] + 24} 452" fill="none" stroke="${C.indigo}" stroke-width="2.2" stroke-opacity=".55" stroke-linecap="round" marker-end="url(#arrowIndigo)"/>`
    }
    nodes += `<path d="M ${x} 322 L ${x} 448" stroke="${C.line}" stroke-width="1.4" stroke-dasharray="4 6"/>`
  })

  const body = `
    ${axis}
    ${nodes}
    ${pill(80, 232, lang === "zh" ? 280 : 300, t.forward, C.teal, C.tealSoft, { h: 32, size: 13 })}
    ${pill(80, 510, lang === "zh" ? 320 : 360, t.backward, C.indigo, C.indigoSoft, { h: 32, size: 13 })}
    ${pill(xs[3] - 78, 388, lang === "zh" ? 220 : 200, t.target, C.rose, C.white, { h: 30, size: 12 })}
    ${pill(xs[5] - 60, 356, lang === "zh" ? 120 : 130, t.known, C.indigo, C.white, { h: 30, size: 12 })}
    ${card(80, 566, 400, 158, {
      title: t.interpTitle,
      body: [t.interp, t.zproc],
      accent: C.indigo,
      fill: C.indigoSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${card(500, 566, 400, 158, {
      title: t.quadTitle,
      body: [t.quad, t.space],
      accent: C.amber,
      fill: C.amberSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${card(920, 566, 420, 158, {
      title: t.orderTitle,
      body: [t.order, t.cost],
      accent: C.rose,
      fill: C.roseSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
  `
  return frame({
    width: 1400,
    height: 770,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function martingaleTraining(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "MARTINGALE LEARNING",
          title: "把方程残差改写成鞅性质",
          subtitle: "残差为零等价于某个过程是鞅；判别网络负责检验，值网络负责修正",
          pdeTitle: "原始表述",
          pdeBody: ["求解算子方程的残差为零", "高维时无法在网格上检验"],
          martTitle: "等价表述",
          martBody: ["构造沿轨道的过程 Mₜ", "残差为零 ⟺ Mₜ 是鞅"],
          testTitle: "判别网络",
          testBody: ["用可测函数族检验鞅性质", "取使违背最大的检验方向"],
          valueTitle: "值网络",
          valueBody: ["最小化被检出的违背量", "无需显式求 inf_u H"],
          controlTitle: "控制网络",
          controlBody: ["在同一循环内更新反馈控制", "不再逐点求解最优化子问题"],
          minmax: "min–max 训练",
          minmaxBody: "值网络与控制网络下降，判别网络上升",
          gain: "为何能进入很高维度",
          gainBody: ["所有量都是沿模拟轨道的期望", "整个算法不出现空间网格"],
        }
      : {
          kicker: "MARTINGALE LEARNING",
          title: "Rewrite the Equation Residual as a Martingale Property",
          subtitle:
            "A vanishing residual is equivalent to a martingale; a test network checks it and a value network fixes it",
          pdeTitle: "Original statement",
          pdeBody: [
            "make the operator residual vanish",
            "unverifiable on a grid in high dimension",
          ],
          martTitle: "Equivalent statement",
          martBody: ["build a path process Mₜ", "residual zero ⟺ Mₜ is a martingale"],
          testTitle: "Test network",
          testBody: ["probes the martingale property", "picks the most violated direction"],
          valueTitle: "Value network",
          valueBody: ["minimises the detected violation", "no explicit inf_u H is needed"],
          controlTitle: "Control network",
          controlBody: ["updates the feedback control", "inside the same training loop"],
          minmax: "min-max training",
          minmaxBody: "value and control descend, the test network ascends",
          gain: "Why high dimensions work",
          gainBody: ["every quantity is an expectation", "along paths, so no grid appears"],
        }

  const body = `
    ${card(40, 210, 380, 170, {
      title: t.pdeTitle,
      body: t.pdeBody,
      accent: C.muted,
      fill: C.slateSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${card(40, 420, 380, 170, {
      title: t.martTitle,
      body: t.martBody,
      accent: C.teal,
      fill: C.tealSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${pathArrow("M 230 386 L 230 414", { color: C.teal, width: 3 })}
    ${section(460, 196, 500, 448, t.minmax, C.indigo, { fill: C.white })}
    ${pill(490, 596, 440, t.minmaxBody, C.indigo, C.indigoSoft, { h: 34, size: 12 })}
    ${card(490, 266, 440, 140, {
      title: t.testTitle,
      body: t.testBody,
      accent: C.rose,
      fill: C.roseSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${card(490, 440, 440, 140, {
      title: t.valueTitle,
      body: t.valueBody,
      accent: C.blue,
      fill: C.blueSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${lineArrow(600, 410, 600, 436, { color: C.rose, width: 2.8 })}
    ${lineArrow(820, 436, 820, 410, { color: C.blue, width: 2.8 })}
    ${card(1000, 266, 360, 140, {
      title: t.controlTitle,
      body: t.controlBody,
      accent: C.amber,
      fill: C.amberSoft,
      align: "center",
      titleSize: 17,
      bodySize: 12,
    })}
    ${lineArrow(968, 336, 992, 336, { color: C.indigo, width: 3 })}
    ${card(1000, 440, 360, 140, {
      title: t.gain,
      body: t.gainBody,
      accent: C.green,
      fill: C.greenSoft,
      align: "center",
      titleSize: 17,
      bodySize: 12,
    })}
    ${lineArrow(968, 510, 992, 510, { color: C.blue, width: 3 })}
  `
  return frame({
    width: 1400,
    height: 700,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function pintDiagonalization(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "DIAGONALIZATION",
          title: "把串行时间递推换成可对角化的时间矩阵",
          subtitle: "代价是特征向量矩阵的条件数：步长比越激进，并行度越好而舍入放大越严重",
          seqTitle: "串行结构",
          seqBody: ["下三角块 Toeplitz 系统", "必须按时间顺序回代"],
          gridTitle: "几何步长",
          gridBody: ["τₙ = τ₁ γⁿ⁻¹，γ > 1", "时间矩阵特征值互不相同"],
          diagTitle: "对角化",
          diagBody: ["B = V D V⁻¹", "D 的元素给出复移位"],
          solveTitle: "并发空间求解",
          solveBody: ["每个移位系统独立", "可放在不同进程上"],
          backTitle: "回到时间层",
          backBody: ["乘回 V 得到全部时间层", "整个过程无外层迭代"],
          tradeTitle: "条件数权衡",
          tradeBody: ["γ 越接近 1，cond(V) 越大", "γ 越大，末端时间步越粗"],
          altTitle: "另一条路线：α-循环预条件",
          altBody: ["等步长下用循环矩阵近似时间矩阵", "对角化由 FFT 完成，条件数受控"],
          use: "两种用法",
          uses: ["直接求解：一次变换即得结果", "预条件：放进 Krylov 迭代做加速"],
        }
      : {
          kicker: "DIAGONALIZATION",
          title: "Replace the Sequential Recurrence by a Diagonalisable Time Matrix",
          subtitle:
            "The price is the eigenvector conditioning: aggressive step growth buys concurrency and amplifies roundoff",
          seqTitle: "Sequential structure",
          seqBody: ["block lower-triangular", "Toeplitz in time"],
          gridTitle: "Geometric steps",
          gridBody: ["τₙ = τ₁ γⁿ⁻¹ with γ > 1", "distinct time eigenvalues"],
          diagTitle: "Diagonalise",
          diagBody: ["B = V D V⁻¹", "D gives complex shifts"],
          solveTitle: "Concurrent solves",
          solveBody: ["each shifted system", "is fully independent"],
          backTitle: "Back to time levels",
          backBody: ["multiply by V", "no outer iteration"],
          tradeTitle: "Conditioning trade-off",
          tradeBody: ["γ near 1 inflates cond(V)", "large γ coarsens the late steps"],
          altTitle: "The other route: α-circulant preconditioning",
          altBody: [
            "uniform steps, a circulant approximation in time",
            "diagonalised by FFT with controlled conditioning",
          ],
          use: "Two ways to use it",
          uses: [
            "direct solver: one transform gives the answer",
            "preconditioner: accelerate a Krylov iteration",
          ],
        }

  const xs = [31, 303, 575, 847, 1119]
  const colors = [C.muted, C.blue, C.indigo, C.teal, C.green]
  const fills = [C.slateSoft, C.blueSoft, C.indigoSoft, C.tealSoft, C.greenSoft]
  const titles = [t.seqTitle, t.gridTitle, t.diagTitle, t.solveTitle, t.backTitle]
  const bodies = [t.seqBody, t.gridBody, t.diagBody, t.solveBody, t.backBody]
  const stages = titles
    .map(
      (title, i) => `${card(xs[i], 214, 250, 160, {
        title,
        body: bodies[i],
        accent: colors[i],
        fill: fills[i],
        step: `${i + 1}`,
        titleSize: 16,
        bodySize: 12,
      })}
      ${i < 4 ? lineArrow(xs[i] + 257, 294, xs[i + 1] - 8, 294, { color: colors[i], width: 3 }) : ""}`,
    )
    .join("")

  const body = `
    ${stages}
    ${card(31, 416, 660, 168, {
      title: t.tradeTitle,
      body: t.tradeBody,
      accent: C.rose,
      fill: C.roseSoft,
      align: "center",
      titleSize: 18,
      bodySize: 14,
    })}
    ${card(709, 416, 660, 168, {
      title: t.altTitle,
      body: t.altBody,
      accent: C.amber,
      fill: C.amberSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${textBlock(700, 626, [t.use], { size: 15, weight: 800, fill: C.ink })}
    ${pill(200, 646, lang === "zh" ? 440 : 470, t.uses[0], C.indigo, C.white, { h: 40, size: 13 })}
    ${pill(lang === "zh" ? 680 : 710, 646, lang === "zh" ? 440 : 450, t.uses[1], C.teal, C.white, { h: 40, size: 13 })}
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

function variableStepEnergy(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "VARIABLE-STEP STABILITY",
          title: "变步长离散的能量论证由二次型正定性支撑",
          subtitle: "步长比进入一个实二次型；该二次型正定，离散能量律与最大值原理才成立",
          problem: "问题",
          problemBody: ["自适应步长破坏常步长的求和技巧", "旧的能量论证不再逐项成立"],
          kernelTitle: "卷积核重排",
          kernelBody: ["用离散正交卷积核抵消历史项", "把多步差分化成单步形式"],
          formTitle: "实二次型",
          formBody: ["历史项汇成 Σ 系数 vⱼ vₖ", "其正定性只依赖步长比序列"],
          energyTitle: "离散能量律",
          energyBody: ["能量单调不增，无需步长上界", "或仅需与界面宽度无关的上界"],
          mbpTitle: "最大值原理",
          mbpBody: ["非线性项的单调性加正定性", "给出逐层的界"],
          ladder: "文献中的三条步长比门槛",
          ladderRows: [
            ["零稳定性", "经典变步长 BDF2 的必要条件", C.blue],
            ["能量稳定性", "离散能量论证给出的更强限制", C.amber],
            ["L² 稳定性", "二次型分析给出的最宽门槛", C.green],
          ],
          note: "同一套工具也用于三阶 BDF、时间滤波 Euler、分数阶 L1 逼近与 IMEX Runge-Kutta",
        }
      : {
          kicker: "VARIABLE-STEP STABILITY",
          title: "Energy Arguments for Variable Steps Rest on a Quadratic Form",
          subtitle:
            "The step ratios enter a real quadratic form; positive definiteness delivers the energy law and the maximum bound",
          problem: "The obstruction",
          problemBody: [
            "adaptive steps break uniform-step telescoping",
            "term-by-term energy arguments fail",
          ],
          kernelTitle: "Convolution kernels",
          kernelBody: [
            "orthogonal kernels cancel the history",
            "a multistep difference becomes one-step",
          ],
          formTitle: "A real quadratic form",
          formBody: [
            "history collects into Σ coefficients vⱼ vₖ",
            "definiteness depends only on step ratios",
          ],
          energyTitle: "Discrete energy law",
          energyBody: [
            "energy decreases with no step-size cap",
            "or a cap independent of interface width",
          ],
          mbpTitle: "Maximum bound principle",
          mbpBody: ["monotone nonlinearity plus definiteness", "gives a bound at every level"],
          ladder: "Three step-ratio thresholds in this literature",
          ladderRows: [
            ["Zero stability", "the classical requirement for variable-step BDF2", C.blue],
            ["Energy stability", "the stronger limit from a discrete energy argument", C.amber],
            ["L² stability", "the widest threshold from the quadratic form", C.green],
          ],
          note: "The same toolkit serves third-order BDF, time-filtered Euler, fractional L1 rules and IMEX Runge-Kutta",
        }

  const rows = t.ladderRows
    .map(
      (row, i) => `${pill(740, 432 + i * 62, 250, row[0], row[2], C.white, { h: 48, size: 14 })}
      ${pill(1006, 432 + i * 62, 360, row[1], row[2], C.white, { h: 48, size: 12 })}`,
    )
    .join("")

  const body = `
    ${card(40, 196, 420, 176, {
      title: t.problem,
      body: t.problemBody,
      accent: C.rose,
      fill: C.roseSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${card(490, 196, 420, 176, {
      title: t.kernelTitle,
      body: t.kernelBody,
      accent: C.indigo,
      fill: C.indigoSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${card(940, 196, 420, 176, {
      title: t.formTitle,
      body: t.formBody,
      accent: C.teal,
      fill: C.tealSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${lineArrow(467, 284, 482, 284, { color: C.rose, width: 3 })}
    ${lineArrow(917, 284, 932, 284, { color: C.indigo, width: 3 })}
    ${textBlock(360, 402, [lang === "zh" ? "由此得到的两条结论" : "The two conclusions it yields"], { size: 15, weight: 800, fill: C.ink })}
    ${textBlock(1053, 402, [t.ladder], { size: 15, weight: 800, fill: C.ink })}
    ${card(40, 424, 330, 176, {
      title: t.energyTitle,
      body: t.energyBody,
      accent: C.green,
      fill: C.greenSoft,
      align: "center",
      titleSize: 17,
      bodySize: 12,
    })}
    ${card(390, 424, 330, 176, {
      title: t.mbpTitle,
      body: t.mbpBody,
      accent: C.blue,
      fill: C.blueSoft,
      align: "center",
      titleSize: 17,
      bodySize: 12,
    })}
    ${rows}
    ${pill(lang === "zh" ? 300 : 240, 638, lang === "zh" ? 800 : 920, t.note, C.muted, C.white, { h: 40, size: 13 })}
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

function unboundedSpectral(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "UNBOUNDED DOMAINS",
          title: "让基函数的衰减与解的衰减一致",
          subtitle: "分数阶算子的解在无穷远只是代数衰减；基函数选错会毁掉谱精度",
          decayTitle: "解的远场行为",
          decayBody: ["局部光滑，远场代数衰减", "衰减指数由算子的阶决定"],
          hermiteTitle: "Hermite 函数",
          hermiteBody: ["自带高斯衰减，适合指数衰减解", "需要引入缩放因子"],
          rationalTitle: "有理基函数",
          rationalBody: ["由代数映射拉回参考区间", "远场按代数速率衰减"],
          symbolTitle: "算子的作用",
          symbolBody: ["分数阶 Laplacian 的符号为 |ξ| 的 2s 次幂", "也可写成超奇异积分"],
          matrixTitle: "离散后的关键量",
          matrixBody: ["基函数在算子下的像", "决定刚度矩阵能否显式算出"],
          matpowTitle: "矩阵的分数次幂",
          matpowBody: ["把 A 的 α 次幂写成有理函数", "或写成围道积分"],
          nodeTitle: "求积节点的代价",
          nodeBody: ["每个节点是一次移位线性求解", "节点之间彼此独立"],
          lesson: "选基的准则",
          lessons: ["匹配远场衰减率", "保持算子像的可积性", "让刚度矩阵结构可利用"],
          algebraic: "代数尾部",
          gaussian: "高斯尾部",
        }
      : {
          kicker: "UNBOUNDED DOMAINS",
          title: "Match the Decay of the Basis to the Decay of the Solution",
          subtitle:
            "Solutions of fractional operators decay only algebraically, and the wrong basis destroys spectral accuracy",
          decayTitle: "Far-field behaviour",
          decayBody: [
            "smooth locally, algebraic at infinity",
            "the rate follows the operator order",
          ],
          hermiteTitle: "Hermite functions",
          hermiteBody: ["Gaussian decay suits exponential tails", "and needs a scaling factor"],
          rationalTitle: "Rational basis",
          rationalBody: ["pulled back by an algebraic map", "so the tail decays algebraically"],
          symbolTitle: "Action of the operator",
          symbolBody: [
            "the fractional Laplacian has symbol |ξ| to the 2s",
            "or a hypersingular integral form",
          ],
          matrixTitle: "The decisive object",
          matrixBody: ["the image of each basis function", "decides if the stiffness is explicit"],
          matpowTitle: "Fractional matrix powers",
          matpowBody: ["write the α power of A as a rational", "function or as a contour integral"],
          nodeTitle: "Cost per quadrature node",
          nodeBody: ["each node is one shifted linear solve", "and the nodes are independent"],
          lesson: "Criteria for choosing a basis",
          lessons: [
            "match the far-field rate",
            "keep the operator image integrable",
            "keep exploitable stiffness structure",
          ],
          algebraic: "algebraic tail",
          gaussian: "Gaussian tail",
        }

  const body = `
    ${section(40, 196, 360, 246, t.decayTitle, C.blue, { fill: C.white })}
    <path d="M 76 410 L 380 410" stroke="${C.line}" stroke-width="1.6"/>
    <path d="M 80 410 L 80 268" stroke="${C.line}" stroke-width="1.6"/>
    <path d="M 88 274 C 168 306, 240 348, 300 372 S 348 390, 374 394" fill="none" stroke="${C.rose}" stroke-width="3.4" stroke-linecap="round"/>
    <path d="M 88 274 C 122 292, 146 384, 196 396 S 300 402, 374 402" fill="none" stroke="${C.blue}" stroke-width="3.4" stroke-linecap="round"/>
    ${pill(86, 416, 132, t.algebraic, C.rose, C.white, { h: 26, size: 11 })}
    ${pill(238, 416, 142, t.gaussian, C.blue, C.white, { h: 26, size: 11 })}
    ${textBlock(220, 480, [t.lesson], { size: 15, weight: 800, fill: C.ink })}
    ${pill(40, 502, 360, t.lessons[0], C.blue, C.white, { h: 38, size: 12 })}
    ${pill(40, 550, 360, t.lessons[1], C.indigo, C.white, { h: 38, size: 12 })}
    ${pill(40, 598, 360, t.lessons[2], C.teal, C.white, { h: 38, size: 12 })}
    ${card(430, 196, 300, 168, {
      title: t.hermiteTitle,
      body: t.hermiteBody,
      accent: C.blue,
      fill: C.blueSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${card(430, 396, 300, 168, {
      title: t.rationalTitle,
      body: t.rationalBody,
      accent: C.rose,
      fill: C.roseSoft,
      align: "center",
      titleSize: 18,
      bodySize: 13,
    })}
    ${card(760, 196, 300, 168, {
      title: t.symbolTitle,
      body: t.symbolBody,
      accent: C.indigo,
      fill: C.indigoSoft,
      align: "center",
      titleSize: 17,
      bodySize: 12,
    })}
    ${card(760, 396, 300, 168, {
      title: t.matrixTitle,
      body: t.matrixBody,
      accent: C.teal,
      fill: C.tealSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${card(1090, 196, 280, 168, {
      title: t.matpowTitle,
      body: t.matpowBody,
      accent: C.amber,
      fill: C.amberSoft,
      align: "center",
      titleSize: 17,
      bodySize: 12,
    })}
    ${card(1090, 396, 280, 168, {
      title: t.nodeTitle,
      body: t.nodeBody,
      accent: C.green,
      fill: C.greenSoft,
      align: "center",
      titleSize: 17,
      bodySize: 12,
    })}
    ${lineArrow(408, 280, 422, 280, { color: C.blue, width: 3 })}
    ${lineArrow(408, 480, 422, 480, { color: C.rose, width: 3 })}
    ${lineArrow(737, 280, 752, 280, { color: C.blue, width: 3 })}
    ${lineArrow(737, 480, 752, 480, { color: C.rose, width: 3 })}
    ${lineArrow(1067, 280, 1082, 280, { color: C.indigo, width: 3 })}
    ${lineArrow(1225, 372, 1225, 388, { color: C.amber, width: 3 })}
  `
  return frame({
    width: 1400,
    height: 660,
    kicker: t.kicker,
    title: t.title,
    subtitle: t.subtitle,
    body,
  })
}

function lowRankDynamics(lang) {
  const t =
    lang === "zh"
      ? {
          kicker: "LOW-RANK DYNAMICS",
          title: "在低秩流形上直接演化随机解",
          subtitle: "基与系数同时随时间变化；规范条件固定了分解的多余自由度",
          ansatzTitle: "分解形式",
          ansatzBody: ["u(t,x,ω) ≈ Σ Uᵢ(t,x) Yᵢ(t,ω)", "秩 R 远小于原始自由度"],
          manifoldTitle: "秩 R 流形",
          manifoldBody: ["所有秩不超过 R 的场构成流形", "真解一般不在流形上"],
          projTitle: "切空间投影",
          projBody: ["把右端投影到切空间", "得到 U 与 Y 的耦合演化方程"],
          gaugeTitle: "规范条件",
          gaugeBody: ["要求基的时间导数与基正交", "消去分解的旋转自由度"],
          errorTitle: "误差来源",
          errorBody: ["最佳秩 R 逼近的距离", "加上切空间投影带来的偏离"],
          riskTitle: "退化风险",
          riskBody: ["最小奇异值趋零时曲率变大", "误差常数随之放大"],
          fixTitle: "两条缓解路线",
          fixList: ["自适应调整秩，及时补充新方向", "对随机方向单独构造演化，减少全局耦合"],
        }
      : {
          kicker: "LOW-RANK DYNAMICS",
          title: "Evolve a Random Solution Directly on a Low-Rank Manifold",
          subtitle:
            "Basis and coefficients both move in time, and a gauge condition removes the redundant freedom",
          ansatzTitle: "The ansatz",
          ansatzBody: ["u(t,x,ω) ≈ Σ Uᵢ(t,x) Yᵢ(t,ω)", "rank R far below the full dimension"],
          manifoldTitle: "The rank-R manifold",
          manifoldBody: [
            "fields of rank at most R form a manifold",
            "the true solution is not on it",
          ],
          projTitle: "Tangent projection",
          projBody: [
            "project the right-hand side onto the tangent",
            "coupled evolution for U and Y",
          ],
          gaugeTitle: "Gauge condition",
          gaugeBody: ["basis time derivatives stay orthogonal", "removing the rotational freedom"],
          errorTitle: "Error sources",
          errorBody: [
            "distance to the best rank-R field",
            "plus the drift from tangent projection",
          ],
          riskTitle: "Degeneracy risk",
          riskBody: [
            "a vanishing smallest singular value",
            "raises curvature and the error constant",
          ],
          fixTitle: "Two mitigations",
          fixList: [
            "adapt the rank and inject new directions",
            "evolve stochastic directions separately",
          ],
        }

  const body = `
    ${card(40, 200, 400, 168, {
      title: t.ansatzTitle,
      body: t.ansatzBody,
      accent: C.blue,
      fill: C.blueSoft,
      align: "center",
      titleSize: 18,
      bodySize: 14,
    })}
    ${section(480, 196, 400, 320, t.manifoldTitle, C.indigo, { fill: C.white })}
    <path d="M 520 430 C 590 340, 700 486, 838 372" fill="none" stroke="${C.indigo}" stroke-width="3.4" stroke-linecap="round"/>
    <circle cx="648" cy="404" r="9" fill="${C.indigo}"/>
    ${textBlock(648, 456, [lang === "zh" ? "当前低秩近似" : "current low-rank state"], { size: 12, weight: 650, fill: C.indigo })}
    <circle cx="700" cy="300" r="9" fill="${C.rose}"/>
    ${textBlock(738, 286, [lang === "zh" ? "真解" : "true solution"], { size: 12, weight: 650, fill: C.rose })}
    ${pathArrow("M 654 396 L 692 308", { color: C.rose, width: 2.4, dashed: true })}
    ${pathArrow("M 648 404 L 792 386", { color: C.teal, width: 2.8 })}
    ${textBlock(770, 356, [lang === "zh" ? "切空间方向" : "tangent direction"], { size: 12, weight: 650, fill: C.teal })}
    ${card(920, 196, 440, 152, {
      title: t.projTitle,
      body: t.projBody,
      accent: C.teal,
      fill: C.tealSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${card(920, 364, 440, 152, {
      title: t.gaugeTitle,
      body: t.gaugeBody,
      accent: C.amber,
      fill: C.amberSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${card(40, 400, 400, 152, {
      title: t.errorTitle,
      body: t.errorBody,
      accent: C.green,
      fill: C.greenSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${card(40, 570, 400, 152, {
      title: t.riskTitle,
      body: t.riskBody,
      accent: C.rose,
      fill: C.roseSoft,
      align: "center",
      titleSize: 17,
      bodySize: 13,
    })}
    ${textBlock(700, 566, [t.fixTitle], { size: 15, weight: 800, fill: C.ink })}
    ${pill(490, 592, lang === "zh" ? 440 : 440, t.fixList[0], C.indigo, C.white, { h: 44, size: 13 })}
    ${pill(950, 592, 410, t.fixList[1], C.teal, C.white, { h: 44, size: 13 })}
    ${lineArrow(447, 284, 472, 284, { color: C.blue, width: 3 })}
    ${lineArrow(888, 272, 912, 272, { color: C.indigo, width: 3 })}
  `
  return frame({
    width: 1400,
    height: 760,
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
  ["tao-zhou-papers", "research-map", researchMap],
  ["tao-zhou-papers", "research-timeline", researchTimeline],
  ["tao-zhou-papers", "sampling-design", samplingDesign],
  ["tao-zhou-papers", "sparse-recovery", sparseRecovery],
  ["tao-zhou-papers", "bayesian-surrogate-loop", bayesianSurrogateLoop],
  ["tao-zhou-papers", "failure-informed-sampling", failureInformedSampling],
  ["tao-zhou-papers", "density-flow-solvers", densityFlowSolvers],
  ["tao-zhou-papers", "operator-learning-uq", operatorLearningUq],
  ["tao-zhou-papers", "fbsde-multistep", fbsdeMultistep],
  ["tao-zhou-papers", "martingale-training", martingaleTraining],
  ["tao-zhou-papers", "pint-diagonalization", pintDiagonalization],
  ["tao-zhou-papers", "variable-step-energy", variableStepEnergy],
  ["tao-zhou-papers", "unbounded-spectral", unboundedSpectral],
  ["tao-zhou-papers", "low-rank-dynamics", lowRankDynamics],
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
